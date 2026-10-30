//! Explicit pairwise-matrix plotting. Labels identify both axes; no inference or conversion.
use crate::{
    cli::Options,
    input::Missing,
    model::{self, Bounds, Content, Figure, Plot},
};
use anyhow::{Context, Result, ensure};
use std::{collections::BTreeSet, io::Read, path::Path};

#[derive(Clone, Debug, PartialEq)]
pub struct Request {
    pub files: Vec<String>,
    pub labels: Vec<String>,
    pub order: Option<String>,
    pub scale: String,
    pub sense: String,
    pub version: u32,
}
impl Default for Request {
    fn default() -> Self {
        Self {
            files: vec![],
            labels: vec![],
            order: None,
            scale: "independent".into(),
            sense: "unspecified".into(),
            version: 1,
        }
    }
}
fn text(path: &str) -> Result<String> {
    ensure!(
        !path.starts_with("1D:") && path != "-" && !path.ends_with(['\'', ']', '}']),
        "RDM inputs must be regular files without AFNI selectors"
    );
    let metadata = std::fs::metadata(path).with_context(|| format!("opening {path}"))?;
    ensure!(
        metadata.is_file() && metadata.len() <= 16 * 1024 * 1024,
        "RDM inputs must be regular files up to 16 MiB"
    );
    let mut text = String::new();
    std::fs::File::open(path)?
        .take(16 * 1024 * 1024 + 1)
        .read_to_string(&mut text)?;
    ensure!(
        text.len() <= 16 * 1024 * 1024,
        "RDM input grew beyond 16 MiB"
    );
    Ok(text)
}
pub fn labels(path: &str) -> Result<Vec<String>> {
    let labels: Vec<String> = text(path)?
        .lines()
        .map(str::trim)
        .filter(|s| !s.is_empty() && !s.starts_with('#'))
        .take(257)
        .map(str::to_owned)
        .collect();
    ensure!(
        (2..=256).contains(&labels.len()),
        "condition labels need 2..256 entries, one label per line"
    );
    ensure!(
        labels
            .iter()
            .all(|s| s.len() <= 256 && !s.chars().any(char::is_control)),
        "condition labels must be at most 256 bytes without tabs/control characters"
    );
    ensure!(
        labels.iter().collect::<BTreeSet<_>>().len() == labels.len(),
        "duplicate condition labels; identities must be unique"
    );
    Ok(labels)
}
pub fn matrix(path: &str, missing: Missing) -> Result<Vec<Vec<f64>>> {
    ensure!(missing != Missing::Mean, "RDMs do not impute values");
    let mut rows: Vec<Vec<f64>> = vec![];
    for line in text(path)?.lines() {
        let line = line.split('#').next().unwrap().trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<_> = line.split_whitespace().take(257).collect();
        ensure!(
            fields.len() <= 256 && rows.len() < 256,
            "RDM matrices support at most 256 rows/columns"
        );
        let row: Vec<_> = fields
            .into_iter()
            .map(|s| {
                let v = if ["na", "n/a", "nan"].contains(&s.to_ascii_lowercase().as_str()) {
                    f64::NAN
                } else {
                    s.parse::<f64>().with_context(|| {
                        format!("RDM matrix must be numeric; invalid cell '{s}'")
                    })?
                };
                ensure!(!v.is_infinite(), "RDM values cannot be infinite");
                ensure!(
                    v.is_finite() || missing == Missing::Gap,
                    "missing RDM value; use -missing gap explicitly"
                );
                Ok(v)
            })
            .collect::<Result<_>>()?;
        rows.push(row);
    }
    let n = rows.len();
    ensure!(
        n >= 2 && rows.iter().all(|r| r.len() == n),
        "RDM matrix must be square with 2..256 conditions; triangle vectors and RSA effect tables are not matrices"
    );
    for (r, row) in rows.iter().enumerate() {
        for (c, other) in rows.iter().enumerate().skip(r + 1) {
            let (a, b) = (row[c], other[r]);
            ensure!(
                (a.is_nan() && b.is_nan())
                    || (a.is_finite()
                        && b.is_finite()
                        && (a - b).abs() <= 1e-6 * a.abs().max(b.abs()).max(1.)),
                "RDM matrix is asymmetric at row {}, column {}; paired missing cells must also match",
                r + 1,
                c + 1
            );
        }
    }
    ensure!(
        rows.iter().flatten().any(|v| v.is_finite()),
        "RDM matrix has no finite values"
    );
    Ok(rows)
}
fn permutation(source: &[String], target: &[String]) -> Result<Vec<usize>> {
    ensure!(
        source.len() == target.len(),
        "condition count mismatch between matrix label sets/order"
    );
    target
        .iter()
        .map(|name| {
            source.iter().position(|s| s == name).with_context(|| {
                format!("condition label mismatch: '{name}' is absent from a matrix's labels")
            })
        })
        .collect()
}
fn range(values: impl Iterator<Item = f64>, zero: bool) -> [f64; 2] {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for v in values.filter(|v| v.is_finite()) {
        lo = lo.min(v);
        hi = hi.max(v);
    }
    if zero {
        let max = lo.abs().max(hi.abs()).max(f64::EPSILON);
        [-max, max]
    } else if lo == hi {
        model::padded(lo, hi)
    } else {
        [lo, hi]
    }
}
pub fn figure(o: &Options) -> Result<Figure> {
    let r = &o.rdm;
    ensure!(r.version == 1, "unsupported RDM preset version");
    ensure!((1..=8).contains(&r.files.len()), "choose 1..8 RDM matrices");
    ensure!(
        r.labels.len() == 1 || r.labels.len() == r.files.len(),
        "supply -rdm-labels FILE once for a common source order, or once per matrix in input order"
    );
    ensure!(
        ["independent", "shared"].contains(&r.scale.as_str()),
        "-rdm-scale must be independent or shared"
    );
    ensure!(
        ["unspecified", "similarity", "dissimilarity"].contains(&r.sense.as_str()),
        "-rdm-kind must be unspecified, similarity, or dissimilarity"
    );
    ensure!(
        ["viridis", "Reds_and_Blues", "Reds_and_Blues_Inv"].contains(&o.cmap.as_str()),
        "unsupported color map"
    );
    let sources: Vec<_> = r.labels.iter().map(|p| labels(p)).collect::<Result<_>>()?;
    let order = r
        .order
        .as_ref()
        .map(|p| labels(p))
        .transpose()?
        .unwrap_or_else(|| sources[0].clone());
    let mut matrices = vec![];
    for (i, file) in r.files.iter().enumerate() {
        let original = matrix(file, o.missing).with_context(|| format!("matrix {file}"))?;
        let source = &sources[if sources.len() == 1 { 0 } else { i }];
        ensure!(
            original.len() == source.len(),
            "{}: {} matrix rows but {} condition labels",
            file,
            original.len(),
            source.len()
        );
        let indices = permutation(source, &order)?;
        matrices.push(
            indices
                .iter()
                .map(|&row| {
                    indices
                        .iter()
                        .map(|&col| original[row][col])
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>(),
        );
    }
    let common = o.clim.or_else(|| {
        (r.scale == "shared")
            .then(|| range(matrices.iter().flatten().flatten().copied(), o.zerocenter))
    });
    let scale_label = if o.clim.is_some() {
        "fixed color scale"
    } else if common.is_some() {
        "shared color scale"
    } else {
        "independent color scale"
    };
    let n = order.len();
    let mut plots = vec![];
    for (i, matrix) in matrices.into_iter().enumerate() {
        let clim = common.unwrap_or_else(|| range(matrix.iter().flatten().copied(), o.zerocenter));
        ensure!(
            clim.iter().all(|v| v.is_finite())
                && clim[0] < clim[1]
                && (clim[1] - clim[0]).is_finite(),
            "RDM color range overflow"
        );
        plots.push(Plot {
            style: Default::default(),
            id: i,
            title: Path::new(&r.files[i])
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .into(),
            note: format!("{n} conditions · {} · {scale_label}", r.sense),
            caption: "".into(),
            xlabel: "Condition".into(),
            ylabel: "Condition".into(),
            content: Content::Heat {
                matrix,
                labels: order.clone(),
                row_labels: order.clone(),
                square_cells: true,
                named_conditions: true,
                cmap: o.cmap.clone(),
                clim,
            },
            bounds: Bounds {
                x: [0., n as f64],
                y: [0., n as f64],
            },
            diagonal: false,
            refs: vec![],
            legend: false,
            qc: Default::default(),
        });
    }
    let mut fig = model::presentation(
        Figure {
            title: o
                .title
                .clone()
                .unwrap_or_else(|| "Representational matrices".into()),
            subtitle: format!(
                "RDM preset v1 · {n} named conditions · {} matrices",
                plots.len()
            ),
            caption: String::new(),
            columns: 2.min(plots.len()),
            plots,
            synthetic: false,
        },
        o,
    )?;
    if fig.caption.is_empty() {
        fig.caption = format!(
            "Conditions aligned by label on both axes · {scale_label}. Original values and diagonal preserved; gray cells indicate missing values. No similarity/distance conversion."
        );
    }
    Ok(fig)
}
pub fn freeze(o: &Options) -> Result<Vec<String>> {
    if o.rdm.files.is_empty() {
        return crate::rsa::freeze(o);
    }
    let mut args = o.arguments.clone();
    for (key, value) in [
        ("-rdm-version", "1"),
        ("-rdm-scale", &o.rdm.scale),
        ("-rdm-kind", &o.rdm.sense),
        ("-cbar", &o.cmap),
    ] {
        if !args.iter().any(|s| s == key) {
            args.extend([key.into(), value.into()]);
        }
    }
    Ok(args)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture(name: &str) -> String {
        format!("{}/examples/rdm/{name}", env!("CARGO_MANIFEST_DIR"))
    }
    fn options(extra: &[&str]) -> Options {
        let mut args = vec![
            "-rdm".into(),
            fixture("model.1D"),
            fixture("neural.1D"),
            "-rdm-labels".into(),
            fixture("conditions.txt"),
            "-rdm-labels".into(),
            fixture("neural-conditions.txt"),
        ];
        args.extend(extra.iter().map(|s| s.to_string()));
        crate::cli::parse(&args).unwrap()
    }
    fn heat(p: &Plot) -> (&Vec<Vec<f64>>, &Vec<String>, [f64; 2]) {
        let Content::Heat {
            matrix,
            labels,
            row_labels,
            clim,
            named_conditions,
            ..
        } = &p.content
        else {
            panic!()
        };
        assert_eq!(labels, row_labels);
        assert!(named_conditions);
        (matrix, labels, *clim)
    }
    #[test]
    fn labels_align_both_axes_shared_scales_and_requested_order_preserve_values() {
        let f = figure(&options(&[
            "-rdm-scale",
            "shared",
            "-rdm-order",
            &fixture("display-order.txt"),
        ]))
        .unwrap();
        let (a, labels, ca) = heat(&f.plots[0]);
        let (b, other, cb) = heat(&f.plots[1]);
        assert_eq!(labels, other);
        assert_eq!(&labels[..2], ["Faces", "Animals"]);
        assert_eq!(a[0][1], 0.18);
        assert_eq!(b[0][1], -0.04);
        assert_eq!(b[1][0], -0.04);
        assert_eq!(ca, [-0.04, 0.95]);
        assert_eq!(ca, cb);
        assert!(
            f.plots[1]
                .hover(1.5, 5.5)
                .contains("Faces × Animals   ·   -0.04000")
        );
        let f = figure(&options(&[])).unwrap();
        assert_eq!(heat(&f.plots[0]).2, [0., 0.92]);
        assert_eq!(heat(&f.plots[1]).2, [-0.04, 0.95]);
        let f = figure(&options(&["-rdm-scale", "shared", "-zerocenter"])).unwrap();
        assert_eq!(heat(&f.plots[0]).2, [-0.95, 0.95]);
        let f = figure(&options(&["-clim", "-1:1", "-rdm-kind", "similarity"])).unwrap();
        assert_eq!(heat(&f.plots[0]).2, [-1., 1.]);
        assert_eq!(heat(&f.plots[1]).0[0][3], -0.04);
        assert!(f.plots[0].note.contains("similarity"));
    }
    #[test]
    fn malformed_labels_and_matrix_structure_fail_without_repair_or_imputation() {
        let d = tempfile::tempdir().unwrap();
        let path = d.path().join("input.1D");
        let file = path.to_str().unwrap();
        for invalid in [
            "0 1\n",
            "0 1\n2 0\n",
            "0 NA\n1 0\n",
            "0 inf\ninf 0\n",
            "NA NA\nNA NA\n",
            "roi A B\nA 0 1\nB 1 0\n",
        ] {
            std::fs::write(&path, invalid).unwrap();
            assert!(matrix(file, Missing::Gap).is_err(), "{invalid}");
        }
        std::fs::write(&path, "1 -0.2 NA\n-0.2 1 0.3\nNA 0.3 1\n").unwrap();
        assert!(matrix(file, Missing::Error).is_err());
        assert!(matrix(file, Missing::Mean).is_err());
        let m = matrix(file, Missing::Gap).unwrap();
        assert_eq!(m[0][0], 1.);
        assert_eq!(m[0][1], -0.2);
        assert!(m[0][2].is_nan());
        for invalid in ["A\nA\n", "A\n", "A\tB\nC\n"] {
            std::fs::write(&path, invalid).unwrap();
            assert!(labels(file).is_err());
        }
        assert!(permutation(&["A".into(), "B".into()], &["A".into(), "C".into()]).is_err());
        std::fs::write(&path, "Faces\nAnimals\n").unwrap();
        assert!(
            figure(&options(&["-rdm-order", file]))
                .unwrap_err()
                .to_string()
                .contains("count mismatch")
        );
        let mut o = options(&[]);
        o.rdm.labels = vec![file.into()];
        assert!(figure(&o).unwrap_err().to_string().contains("matrix rows"));
        o.rdm.labels = vec![];
        assert!(figure(&o).is_err());
        let mut o = options(&[]);
        o.rdm.scale = "maybe".into();
        assert!(figure(&o).is_err());
        let mut o = options(&[]);
        o.rdm.sense = "distance-conversion".into();
        assert!(figure(&o).is_err());
    }
    #[test]
    fn matrix_size_tolerance_and_label_uniqueness_are_bounded() {
        let d = tempfile::tempdir().unwrap();
        let path = d.path().join("matrix.1D");
        let file = path.to_str().unwrap();
        std::fs::write(&path, "0 0.1\n0.10000001 0\n").unwrap();
        let m = matrix(file, Missing::Error).unwrap();
        assert_ne!(m[0][1], m[1][0]);
        std::fs::write(&path, (0..257).map(|_| "0").collect::<Vec<_>>().join(" ")).unwrap();
        assert!(matrix(file, Missing::Error).is_err());
        std::fs::write(&path, "# conditions\nFace left\nFace right\n").unwrap();
        assert_eq!(labels(file).unwrap(), ["Face left", "Face right"]);
        std::fs::write(&path, format!("{}\nB\n", "a".repeat(257))).unwrap();
        assert!(labels(file).is_err());
    }
    #[test]
    fn project_tracks_matrix_labels_and_order_and_replays_relocated_pixels() {
        let d = tempfile::tempdir().unwrap();
        for name in [
            "model.1D",
            "neural.1D",
            "conditions.txt",
            "neural-conditions.txt",
            "display-order.txt",
        ] {
            std::fs::copy(fixture(name), d.path().join(name)).unwrap();
        }
        let mut o = options(&[
            "-rdm-scale",
            "shared",
            "-rdm-order",
            &fixture("display-order.txt"),
        ]);
        for i in o.source_args.clone() {
            o.arguments[i] = d
                .path()
                .join(Path::new(&o.arguments[i]).file_name().unwrap())
                .to_str()
                .unwrap()
                .into();
        }
        o = crate::cli::parse(&o.arguments).unwrap();
        assert_eq!(o.source_args.len(), 5);
        let f = figure(&o).unwrap();
        let origin = crate::document::Origin::capture(&o).unwrap();
        let state = crate::document::Snapshot {
            presentation: crate::document::Presentation::capture(&f),
            views: f.plots.iter().map(|p| p.bounds).collect(),
            width: 1500,
            height: 1000,
            full_plot: false,
        };
        let recipe = d.path().join("figure.json");
        crate::document::save(&recipe, &origin, &state, false).unwrap();
        let json = std::fs::read_to_string(&recipe).unwrap();
        assert!(json.contains("-rdm-version"));
        let moved = d.path().join("moved");
        std::fs::create_dir(&moved).unwrap();
        for name in [
            "model.1D",
            "neural.1D",
            "conditions.txt",
            "neural-conditions.txt",
            "display-order.txt",
            "figure.json",
        ] {
            std::fs::rename(d.path().join(name), moved.join(name)).unwrap();
        }
        let loaded = crate::document::load(&moved.join("figure.json")).unwrap();
        let renderer = crate::render::Renderer::new();
        let pixels = |f: &Figure| {
            renderer
                .raster(
                    &renderer
                        .scene(f, &state.views, state.width, state.height)
                        .unwrap(),
                )
                .unwrap()
        };
        assert_eq!(pixels(&f).data(), pixels(&loaded.figure).data());
        std::fs::write(moved.join("conditions.txt"), "changed\n").unwrap();
        assert!(crate::document::load(&moved.join("figure.json")).is_err());
    }
}
