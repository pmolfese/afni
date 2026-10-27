//! Independent x-range annotations and AFNI design-matrix activity (not event timing).
use crate::{
    input,
    model::{COLORS, color},
};
use anyhow::{Context, Result, ensure};
use std::collections::{BTreeMap, HashMap};

#[derive(Clone, Debug, PartialEq)]
pub struct Shade {
    pub x: [f64; 2],
    pub color: String,
    pub alpha: f64,
    pub label: Option<String>,
}
impl Shade {
    pub fn new(x: [f64; 2], index: usize) -> Self {
        Self {
            x,
            color: COLORS[index % COLORS.len()].into(),
            alpha: 0.18,
            label: None,
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct Design {
    pub file: Option<String>,
    pub colors: Vec<String>,
    pub threshold: f64,
}
impl Default for Design {
    fn default() -> Self {
        Self {
            file: None,
            colors: vec![],
            threshold: 0.05,
        }
    }
}
fn integers(s: &str) -> Result<Vec<i64>> {
    let mut out = vec![];
    for v in s
        .split(|c: char| c == ',' || c.is_whitespace())
        .filter(|s| !s.is_empty())
    {
        if let Some((a, b)) = v.split_once("..") {
            let (a, b) = (a.parse::<i64>()?, b.parse::<i64>()?);
            ensure!(
                b >= a && b.checked_sub(a).is_some_and(|n| n < 1_000_000),
                "invalid or excessive X-matrix integer range"
            );
            out.extend(a..=b);
        } else if let Some((n, value)) = v.split_once('@') {
            let n = n.parse::<usize>()?;
            ensure!(n > 0 && n <= 1_000_000, "invalid X-matrix repetition count");
            out.extend(std::iter::repeat_n(value.parse::<i64>()?, n));
        } else {
            out.push(v.parse()?);
        }
        ensure!(out.len() <= 1_000_000, "excessive X-matrix metadata length");
    }
    Ok(out)
}
fn metadata(text: &str) -> Result<HashMap<String, String>> {
    let mut out = HashMap::new();
    for line in text.lines().filter_map(|l| l.trim().strip_prefix('#')) {
        if let Some((key, value)) = line.trim().split_once('=') {
            let key = key.trim();
            if [
                "ColumnLabels",
                "ColumnGroups",
                "StimLabels",
                "StimBots",
                "StimTops",
                "GoodList",
                "NRowFull",
            ]
            .contains(&key)
            {
                let value = value.trim();
                let value = value
                    .strip_prefix('"')
                    .and_then(|s| s.strip_suffix('"'))
                    .context("X-matrix metadata values must be quoted")?;
                ensure!(
                    out.insert(key.into(), value.into()).is_none(),
                    "duplicate X-matrix metadata {key}"
                );
            }
        }
    }
    Ok(out)
}
pub fn design_spans(o: &Design, x: &[f64]) -> Result<Vec<Shade>> {
    let Some(file) = &o.file else {
        return Ok(vec![]);
    };
    ensure!(
        !file.ends_with([']', '}', '\'']) && file != "-" && !file.starts_with("1D:"),
        "shade-xmat needs an unselected X.mat.1D file with its metadata"
    );
    let text =
        std::fs::read_to_string(file).with_context(|| format!("reading shading matrix {file}"))?;
    let meta = metadata(&text)?;
    let table = input::read(file)?;
    let nc = table[0].len();
    let groups = integers(meta.get("ColumnGroups").context(
        "shade-xmat requires ColumnGroups to distinguish conditions from baseline/motion",
    )?)?;
    ensure!(
        groups.len() == nc,
        "X-matrix ColumnGroups count must match data columns"
    );
    let labels: Vec<_> = meta
        .get("ColumnLabels")
        .context("shade-xmat requires ColumnLabels")?
        .split(';')
        .map(str::trim)
        .collect();
    ensure!(
        labels.len() == nc && labels.iter().all(|s| !s.is_empty()),
        "X-matrix ColumnLabels count must match data columns"
    );
    let full = meta
        .get("NRowFull")
        .map(|s| s.parse::<usize>())
        .transpose()?
        .unwrap_or(table.len());
    ensure!(
        full > 0 && x.len() == full,
        "shade-xmat requires the full {full}-sample plotted timeline; use full data and X.nocensor.xmat.1D when available"
    );
    let good = if let Some(g) = meta.get("GoodList") {
        let values = integers(g)?;
        ensure!(
            values.iter().all(|&v| v >= 0 && (v as u64) < full as u64),
            "X-matrix GoodList index outside full timeline"
        );
        values.into_iter().map(|v| v as usize).collect::<Vec<_>>()
    } else {
        ensure!(
            table.len() == full,
            "censored X-matrix needs GoodList to locate rows"
        );
        (0..full).collect()
    };
    // Uncensored AFNI matrices may retain GoodList metadata but contain every row.
    ensure!(
        good.windows(2).all(|w| w[0] < w[1]),
        "X-matrix GoodList must be increasing and unique"
    );
    let row_ids = if table.len() == full {
        (0..full).collect::<Vec<_>>()
    } else {
        ensure!(
            good.len() == table.len(),
            "X-matrix GoodList count must match matrix rows"
        );
        good
    };
    ensure!(
        x.iter().all(|v| v.is_finite()) && x.windows(2).all(|w| w[1] > w[0]),
        "shade-xmat requires increasing finite x coordinates"
    );
    let mut edges = vec![if full == 1 {
        x[0] - 0.5
    } else {
        x[0] - (x[1] - x[0]) / 2.
    }];
    edges.extend(x.windows(2).map(|w| w[0] + (w[1] - w[0]) / 2.));
    edges.push(if full == 1 {
        x[0] + 0.5
    } else {
        x[full - 1] + (x[full - 1] - x[full - 2]) / 2.
    });
    ensure!(
        edges.iter().all(|v| v.is_finite()),
        "X-matrix shading boundaries overflow"
    );
    let mut columns: BTreeMap<i64, Vec<usize>> = BTreeMap::new();
    for (c, &g) in groups.iter().enumerate() {
        if g > 0 {
            columns.entry(g).or_default().push(c);
        }
    }
    ensure!(
        !columns.is_empty(),
        "X-matrix contains no positive condition groups"
    );
    let palette = o
        .colors
        .iter()
        .map(|s| color(s))
        .collect::<Result<Vec<_>>>()?;
    let mut names: HashMap<i64, String> = HashMap::new();
    if let (Some(names_text), Some(bots), Some(tops)) = (
        meta.get("StimLabels"),
        meta.get("StimBots"),
        meta.get("StimTops"),
    ) {
        let ns: Vec<_> = names_text.split(';').map(str::trim).collect();
        let bots = integers(bots)?;
        let tops = integers(tops)?;
        ensure!(
            ns.len() == bots.len() && ns.len() == tops.len(),
            "inconsistent X-matrix stimulus metadata"
        );
        for ((name, b), t) in ns.into_iter().zip(bots).zip(tops) {
            ensure!(
                b >= 0 && t >= b && (t as usize) < nc,
                "X-matrix stimulus column range is invalid"
            );
            if groups[b as usize] > 0 {
                names
                    .entry(groups[b as usize])
                    .or_insert_with(|| name.into());
            }
        }
    }
    let mut spans = vec![];
    for (index, (&group, cols)) in columns.iter().enumerate() {
        let fallback = labels[cols[0]].split('#').next().unwrap_or(labels[cols[0]]);
        let name = names
            .get(&group)
            .cloned()
            .unwrap_or_else(|| fallback.into());
        let label = format!("{name} (model)");
        let peak = table
            .iter()
            .flat_map(|r| cols.iter().map(move |&c| r[c].abs()))
            .fold(0f64, f64::max);
        if peak == 0. {
            continue;
        }
        let mut active = vec![false; full];
        for (r, &sample) in table.iter().zip(&row_ids) {
            active[sample] = cols.iter().any(|&c| r[c].abs() > peak * o.threshold);
        }
        let c = if palette.is_empty() {
            COLORS[index % COLORS.len()].into()
        } else {
            palette[index % palette.len()].clone()
        };
        let mut i = 0;
        while i < full {
            if !active[i] {
                i += 1;
                continue;
            }
            let start = i;
            while i < full && active[i] {
                i += 1;
            }
            spans.push(Shade {
                x: [edges[start], edges[i]],
                color: c.clone(),
                alpha: 0.18,
                label: Some(label.clone()),
            });
        }
    }
    eprintln!(
        "++ X-matrix shading: modeled response > {:.1}% of each condition's peak magnitude; not stimulus on/off timing. {} missing matrix rows remain unshaded.",
        o.threshold * 100.,
        full - table.len()
    );
    Ok(spans)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn design_file(d: &std::path::Path, text: &str) -> Design {
        let p = d.join("X.mat.1D");
        std::fs::write(&p, text).unwrap();
        Design {
            file: Some(p.to_str().unwrap().into()),
            ..Default::default()
        }
    }
    #[test]
    fn afni_metadata_runs_and_repeats() {
        assert_eq!(
            integers("2@-1,0,2@1,2..4").unwrap(),
            [-1, -1, 0, 1, 1, 2, 3, 4]
        );
        assert!(integers("0@1").is_err());
        assert!(integers("999999999@1").is_err());
    }
    #[test]
    fn condition_groups_ignore_nuisance_union_bases_and_preserve_missing_rows() {
        let d = tempfile::tempdir().unwrap();
        let mut o = design_file(
            d.path(),
            "# ColumnLabels = \"base;motion;A#0;A#1;B#0\"\n# ColumnGroups = \"-1,0,2@1,2\"\n# NRowFull = \"5\"\n# GoodList = \"0,1,3,4\"\n# StimLabels = \"Alpha;Beta\"\n# StimBots = \"2,4\"\n# StimTops = \"3,4\"\n1 9 0 0 0\n1 9 2 0 0\n1 9 0 -2 4\n1 9 0 0 4\n",
        );
        o.colors = vec!["blue".into(), "orange".into()];
        let spans = design_spans(&o, &[0., 2., 4., 6., 8.]).unwrap();
        assert_eq!(spans.len(), 3);
        assert_eq!(spans[0].x, [1., 3.]);
        assert_eq!(spans[1].x, [5., 7.]);
        assert_eq!(spans[0].label.as_deref(), Some("Alpha (model)"));
        assert_eq!(spans[2].x, [5., 9.]);
        assert_ne!(spans[0].color, spans[2].color);
        assert!(design_spans(&o, &[0., 1., 2., 3.]).is_err());
    }
    #[test]
    fn malformed_or_ambiguous_matrix_metadata_fails() {
        let d = tempfile::tempdir().unwrap();
        for text in [
            "# ColumnLabels = \"A\"\n1\n2\n",
            "# ColumnLabels = \"A\"\n# ColumnGroups = \"1,2\"\n1\n2\n",
            "# ColumnLabels = \"A\"\n# ColumnGroups = \"1\"\n# GoodList = \"1,1\"\n1\n2\n",
            "# ColumnLabels = \"A\"\n# ColumnGroups = \"0\"\n1\n2\n",
        ] {
            let o = design_file(d.path(), text);
            assert!(design_spans(&o, &[0., 1.]).is_err());
        }
    }
    #[test]
    fn real_deconvolve_fixture_and_adjustable_threshold() {
        let o = Design {
            file: Some(format!(
                "{}/examples/X.demo.mat.1D",
                env!("CARGO_MANIFEST_DIR")
            )),
            ..Default::default()
        };
        let x = (0..120).map(|i| i as f64).collect::<Vec<_>>();
        let spans = design_spans(&o, &x).unwrap();
        assert_eq!(spans.len(), 4);
        assert_eq!(spans[0].label.as_deref(), Some("Task_A (model)"));
        let high = design_spans(
            &Design {
                threshold: 0.8,
                ..o
            },
            &x,
        )
        .unwrap();
        assert!(
            high.iter().map(|s| s.x[1] - s.x[0]).sum::<f64>()
                < spans.iter().map(|s| s.x[1] - s.x[0]).sum::<f64>()
        );
    }
}
