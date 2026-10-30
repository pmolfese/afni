//! Named fitted/commonality quantities from 3dRSA's ordinary table writer.
//! Projection reuses plotting geometry, never the scientific meaning of ordinary beta/r.
use crate::{
    cli::Options,
    model::Figure,
    rsa::{self, Table},
};
use anyhow::{Result, ensure};
use std::collections::BTreeMap;
pub const FAMILIES: [(&str, &str); 6] = [
    ("effects", "Primary effects"),
    ("fit", "Cross-validated accuracy (cvR)"),
    ("weights", "Mean normalized fold weights"),
    ("fit-contrast", "Held-fold Fisher-z accuracy difference"),
    ("commonality", "Raw commonality (R² decomposition)"),
    ("partial-r2", "Conditional partial R²"),
];
pub fn component(t: &Table, name: &str) -> bool {
    ["uniq_", "common_", "partialR2_"]
        .iter()
        .any(|p| name.starts_with(p))
        && t.columns.contains(&format!("{name}_p"))
        && t.columns.contains(&format!("{name}_q"))
}
pub fn family(t: &Table, index: usize) -> Option<&'static str> {
    if t.long.is_some() {
        return None;
    }
    let name = t.columns.get(index)?;
    if name.ends_with("_cvR") {
        Some("fit")
    } else if name.ends_with("_cvDiff") {
        Some("fit-contrast")
    } else if t
        .columns
        .iter()
        .filter_map(|s| s.strip_suffix("_cvR"))
        .any(|m| name.starts_with(&format!("{m}_w_")))
    {
        Some("weights")
    } else if component(t, name) {
        Some(if name.starts_with("partialR2_") {
            "partial-r2"
        } else {
            "commonality"
        })
    } else {
        None
    }
}
pub fn validate(t: &Table) -> Result<()> {
    for (i, name) in t.columns.iter().enumerate() {
        match family(t, i) {
            Some("fit" | "fit-contrast") => {
                let base = name
                    .strip_suffix("_cvR")
                    .or_else(|| name.strip_suffix("_cvDiff"))
                    .unwrap();
                ensure!(
                    t.columns.contains(&format!("{base}_cvP"))
                        && t.columns.contains(&format!("{base}_cvQ")),
                    "fitted quantity {name} lacks cvP/cvQ"
                );
            }
            Some("weights") => ensure!(!name.ends_with("_w_"), "empty fitted component name"),
            _ => {}
        }
    }
    Ok(())
}
fn project(o: &Options, t: &Table) -> Result<(Options, Table, Vec<String>, &'static str)> {
    let label = FAMILIES
        .iter()
        .find(|(f, _)| *f == o.rsa.family && *f != "effects")
        .map(|(_, l)| *l);
    ensure!(
        label.is_some(),
        "unknown RSA family; choose effects, fit, weights, fit-contrast, commonality, partial-r2"
    );
    ensure!(
        t.long.is_none(),
        "fit/commonality families require a wide 3dRSA table"
    );
    ensure!(
        o.rsa.models.is_none() && o.rsa.pair.is_none() && o.rsa.effect == "auto",
        "use -rsa-measures for this family; -rsa-models/-rsa-effect/-rsa-pair select primary effects"
    );
    let available: Vec<_> = t
        .columns
        .iter()
        .enumerate()
        .filter(|(i, _)| family(t, *i) == Some(o.rsa.family.as_str()))
        .map(|(_, s)| s.clone())
        .collect();
    let chosen = o.rsa.measures.clone().unwrap_or_else(|| available.clone());
    ensure!(
        !chosen.is_empty() && chosen.len() <= 8 && chosen.iter().all(|s| available.contains(s)),
        "choose 1..8 {} measures with -rsa-measures (available: {})",
        o.rsa.family,
        available.join(",")
    );
    let mut proxy = o.clone();
    proxy.rsa.models = Some(chosen.clone());
    proxy.rsa.effect = "b".into();
    proxy.rsa.family = "effects".into();
    proxy.rsa.measures = None;
    if proxy.rsa.plots.is_none() {
        proxy.rsa.plots = Some(vec!["heatmap".into()]);
    }
    let mut table = t.clone();
    table.models = chosen.clone();
    table.effects = BTreeMap::new();
    for name in &chosen {
        let i = t.columns.iter().position(|s| s == name).unwrap();
        table.effects.insert((name.clone(), "b".into()), i);
    }
    // Supplied commonality CIs already have quantity_bootLo/Hi names. Fit/weights do not.
    if ["fit", "weights", "fit-contrast"].contains(&o.rsa.family.as_str()) {
        ensure!(
            !proxy
                .rsa
                .plots
                .as_ref()
                .unwrap()
                .iter()
                .any(|m| m == "intervals")
                && o.rsa.ci == "auto",
            "3dRSA does not supply CI endpoints for fitted accuracy, fit contrasts or mean weights; uncertainty is not inferred"
        );
    }
    Ok((proxy, table, chosen, label.unwrap()))
}
pub fn freeze(o: &Options, t: &Table) -> Result<Vec<String>> {
    let (proxy, table, chosen, _) = project(o, t)?;
    let resolved = rsa::freeze_table(&proxy, &table)?;
    let mut args = o.arguments.clone();
    let mut flags = vec![
        ("-rsa-version", "1".into()),
        ("-rsa-family", o.rsa.family.clone()),
        ("-rsa-measures", chosen.join(",")),
    ];
    for flag in ["-rsa-rois", "-rsa-plots", "-rsa-ci"] {
        if let Some(i) = resolved.iter().position(|s| s == flag) {
            flags.push((flag, resolved[i + 1].clone()));
        }
    }
    for (flag, value) in flags {
        if let Some(i) = args.iter().position(|s| s == flag) {
            args[i + 1] = value;
        } else {
            args.extend([flag.into(), value]);
        }
    }
    Ok(args)
}
pub fn figure(o: &Options, t: &Table) -> Result<Figure> {
    let (proxy, table, _, label) = project(o, t)?;
    let mut figure = rsa::figure_table(&proxy, &table, Some(label))?;
    if o.caption.is_empty() {
        figure.caption=match o.rsa.family.as_str(){
            "weights"=>"Supplied mean fold weights, normalized by 3dRSA; descriptive component allocation, not inferential evidence. No rescaling or intervals inferred.",
            "commonality"=>"Supplied raw R² decomposition components. Negative shared components indicate suppression and are retained. Conditional partial R² is a separate family; no clipping, normalization, sums or significance filtering.",
            "partial-r2"=>"Supplied conditional partial R²: unique contribution relative to variance left by competitors. These quantities are not additive raw decomposition regions. No clipping, normalization or inference recomputation.",
            "fit-contrast"=>"Supplied paired held-fold Fisher-z accuracy differences. No conversion, inference recomputation or significance filtering.",
            _=>"Supplied cross-validated accuracy (cvR), on the source scale. Fold construction is determined by 3dRSA; no refitting, conversion or significance filtering.",
        }.into();
        for meta in t.metadata.iter().filter(|m| {
            m.starts_with("fitted") || m.starts_with("commonality") || m.starts_with("metric")
        }) {
            figure.caption.push_str(&format!(" {meta}."));
        }
    }
    figure.subtitle = figure.subtitle.replace(" models", " quantities");
    for plot in &mut figure.plots {
        if matches!(plot.content, crate::model::Content::Heat { .. }) {
            plot.xlabel = "Quantity".into();
            plot.note = plot.note.replace(" models", " quantities");
            if !o.panel_titles.iter().any(|(i, _)| *i == plot.id + 1) {
                plot.title = format!("ROI × quantity · {label}");
            }
        }
    }
    // Preserve selected effect semantics in all labels, including user-customized titles.
    if o.title.is_none() {
        figure.title = format!("RSA · {label}");
    }
    Ok(figure)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> String {
        format!(
            "{}/examples/rsa/fitted-commonality.rsa.1D",
            env!("CARGO_MANIFEST_DIR")
        )
    }
    fn options(family: &str, flags: &[&str]) -> Options {
        let mut args = vec![
            "-rsa".into(),
            fixture(),
            "-rsa-family".into(),
            family.into(),
        ];
        args.extend(flags.iter().map(|s| s.to_string()));
        crate::cli::parse(&args).unwrap()
    }
    #[test]
    fn families_keep_scales_negative_regions_and_descriptive_weights_separate() {
        for (family, label, n) in [
            ("fit", "Cross-validated", 2),
            ("weights", "fold weights", 2),
            ("fit-contrast", "Fisher-z", 1),
            ("commonality", "decomposition", 3),
            ("partial-r2", "Conditional", 2),
        ] {
            let o = options(family, &["-rsa-rois", "12,3"]);
            let f = rsa::figure(&o).unwrap();
            assert!(f.title.contains(label));
            let crate::model::Content::Heat { matrix, labels, .. } = &f.plots[0].content else {
                panic!()
            };
            assert_eq!(matrix.len(), 2);
            assert_eq!(labels.len(), n);
            if family == "commonality" {
                assert_eq!(matrix[0][2], -0.06);
                assert!(f.caption.contains("suppression"));
            }
            if family == "weights" {
                assert_eq!(matrix[0], [0.65, 0.35]);
                assert!(f.caption.contains("descriptive"));
            }
            let frozen = rsa::freeze(&o).unwrap();
            assert!(frozen.iter().any(|s| s == "-rsa-measures"));
            assert!(!frozen.iter().any(|s| s == "-rsa-models"));
            let replay = rsa::figure(&crate::cli::parse(&frozen).unwrap()).unwrap();
            assert_eq!(
                format!("{:?}", f.plots[0].content),
                format!("{:?}", replay.plots[0].content)
            );
        }
        assert!(rsa::figure(&options("fit", &["-rsa-measures", "combined_cvP"])).is_err());
        assert!(rsa::figure(&options("commonality", &["-rsa-measures", "partialR2_A"])).is_err());
        assert!(rsa::figure(&options("weights", &["-rsa-plots", "intervals"])).is_err());
        assert!(rsa::figure(&options("fit", &["-rsa-effect", "r"])).is_err());
    }
    #[test]
    fn commonality_supplied_bounds_and_triple_labels_remain_exact() {
        let f = rsa::figure(&options(
            "commonality",
            &[
                "-rsa-measures",
                "common_A_B",
                "-rsa-plots",
                "intervals",
                "-rsa-rois",
                "12",
            ],
        ))
        .unwrap();
        let crate::model::Content::Intervals { rows, method, .. } = &f.plots[0].content else {
            panic!()
        };
        assert_eq!(rows[0].values, Some([-0.06, -0.10, 0.01]));
        assert!(method.contains("90%"));
        let mut o = options("commonality", &[]);
        o.rsa.file = Some(format!(
            "{}/examples/rsa/commonality-triple.rsa.1D",
            env!("CARGO_MANIFEST_DIR")
        ));
        let f = rsa::figure(&o).unwrap();
        let crate::model::Content::Heat { labels, matrix, .. } = &f.plots[0].content else {
            panic!()
        };
        assert_eq!(labels.len(), 7);
        assert!(labels.contains(&"common_A_B_not_C".into()));
        assert!(matrix[0][6] < 0.);
        let source = std::fs::read_to_string(fixture()).unwrap();
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("bad.1D");
        std::fs::write(&p, source.replace("combined_cvQ", "unknown")).unwrap();
        assert!(rsa::read(&p).is_err());
    }
}
