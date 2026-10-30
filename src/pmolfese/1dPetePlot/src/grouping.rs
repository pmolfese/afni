//! Explicit categorical series partitioning. Never aggregate, interpolate, or join rows.
use crate::{
    input,
    model::{COLORS, Series},
};
use anyhow::{Context, Result, ensure};
use std::io::Read;

pub fn split(series: Vec<Series>, file: &str, column: &str) -> Result<Vec<Series>> {
    ensure!(
        !file.contains('{') && !file.ends_with('\''),
        "group-by does not support row selection or transpose; provide a table in the intended row order"
    );
    let name = if file.ends_with(']') {
        &file[..file.rfind('[').context("invalid column selector")?]
    } else {
        file
    };
    ensure!(
        name != "-" && !name.starts_with("1D:"),
        "group-by requires a regular .1D or headered TSV file"
    );
    ensure!(
        std::fs::metadata(name)?.len() <= 16 * 1024 * 1024,
        "group-by supports tables up to 16 MiB"
    );
    let mut text = String::new();
    std::fs::File::open(name)?
        .take(16 * 1024 * 1024 + 1)
        .read_to_string(&mut text)?;
    ensure!(
        text.len() <= 16 * 1024 * 1024,
        "group table grew beyond 16 MiB"
    );
    let tsv = name.to_ascii_lowercase().ends_with(".tsv");
    let mut header = vec![];
    let mut labels: Vec<String> = vec![];
    let mut identities = vec![];
    let mut width = None;
    for line in text.lines() {
        let line = line.split('#').next().unwrap_or("");
        if line.trim().is_empty() {
            continue;
        }
        let words: Vec<_> = if tsv {
            line.split('\t').map(str::trim).collect()
        } else {
            line.split_whitespace().collect()
        };
        if tsv && header.is_empty() {
            header = words.clone();
            width = Some(words.len());
            continue;
        }
        ensure!(width.is_none_or(|n| words.len() == n), "ragged group table");
        width = Some(words.len());
        let c = if let Some(i) = header.iter().position(|h| *h == column) {
            i
        } else {
            let selected = input::indices(column, words.len())?;
            ensure!(
                selected.len() == 1,
                "group-by needs exactly one original column"
            );
            selected[0]
        };
        let value = words[c];
        ensure!(
            !["", "na", "n/a", "nan", "inf", "-inf"].contains(&value.to_ascii_lowercase().as_str())
                && value.len() <= 512
                && !value.parse::<f64>().is_ok_and(|v| !v.is_finite()),
            "group labels must be nonmissing and at most 512 bytes; grouping never imputes labels"
        );
        let group = if let Some(i) = labels.iter().position(|s| s == value) {
            i
        } else {
            ensure!(labels.len() < 16, "group-by supports at most 16 groups");
            labels.push(value.into());
            labels.len() - 1
        };
        identities.push(group);
    }
    ensure!(
        !labels.is_empty() && series.iter().all(|s| s.points.len() == identities.len()),
        "group rows must match every plotted sample"
    );
    ensure!(
        series.len() * labels.len() <= 16,
        "group-by supports at most 16 total series (Y columns × groups)"
    );
    let mut output = vec![];
    for source in series {
        for (group, label) in labels.iter().enumerate() {
            let mut s = source.clone();
            s.label = format!("{} · {label}", source.label);
            s.color = COLORS[group % COLORS.len()].into();
            for (row, &identity) in identities.iter().enumerate() {
                if identity != group {
                    s.points[row][1] = f64::NAN;
                    if let Some(band) = &mut s.band {
                        band[row] = [f64::NAN; 2];
                    }
                }
            }
            ensure!(
                s.points.iter().any(|p| p[1].is_finite()),
                "group '{}' has no finite observations for '{}'",
                label,
                source.label
            );
            output.push(s);
        }
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use crate::{
        cli,
        model::{self, Content},
    };
    #[test]
    fn groups_preserve_categories_rows_and_share_histogram_edges() {
        let d = tempfile::tempdir().unwrap();
        let path = d.path().join("table.tsv");
        std::fs::write(&path, "x\ty\tgroup\n0\t1\tB\n1\t2\tA\n2\t4\tB\n3\t5\tA\n").unwrap();
        for kind in ["-histogram", "-scatter"] {
            let mut args = vec![
                "-infiles".into(),
                format!("{}[y]", path.display()),
                "-one".into(),
                "-group-by".into(),
                "group".into(),
                kind.into(),
            ];
            if kind == "-scatter" {
                args.extend(["-xfile".into(), format!("{}[x]", path.display())]);
            }
            let fig = model::from_options(&cli::parse(&args).unwrap()).unwrap();
            match &fig.plots[0].content {
                Content::Histogram { series, edges, .. } => {
                    assert_eq!(edges[0], edges[1]);
                    assert_eq!(series[0].points.iter().map(|p| p[1]).sum::<f64>(), 2.);
                    assert_eq!(series[1].points.iter().map(|p| p[1]).sum::<f64>(), 2.);
                    assert!(series[0].label.ends_with(" · B"));
                }
                Content::Scatter(s) => {
                    assert_eq!(s[0].points[2], [2., 4.]);
                    assert!(s[0].points[1][1].is_nan());
                    assert_eq!(s[1].points[1], [1., 2.]);
                }
                _ => panic!(),
            }
        }
        for values in ["0\t1\tNA\n1\t2\tA", "0\t1\t1e999\n1\t2\tA"] {
            std::fs::write(&path, format!("x\ty\tgroup\n{values}\n")).unwrap();
            let o = cli::parse(&[
                "-infiles".into(),
                format!("{}[y]", path.display()),
                "-one".into(),
                "-group-by".into(),
                "group".into(),
            ])
            .unwrap();
            assert!(model::from_options(&o).is_err());
        }
        std::fs::write(
            &path,
            format!(
                "x\ty\tgroup\n{}",
                (0..17)
                    .map(|i| format!("{i}\t1\tG{i}\n"))
                    .collect::<String>()
            ),
        )
        .unwrap();
        let o = cli::parse(&[
            "-infiles".into(),
            format!("{}[y]", path.display()),
            "-one".into(),
            "-group-by".into(),
            "group".into(),
        ])
        .unwrap();
        assert!(model::from_options(&o).is_err());
    }
}
