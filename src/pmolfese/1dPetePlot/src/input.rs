use anyhow::{Context, Result, bail, ensure};
use std::io::Read;

pub fn indices(spec: &str, n: usize) -> Result<Vec<usize>> {
    ensure!(n > 0, "cannot select from an empty axis");
    let index = |v: &str| -> Result<usize> {
        let i = if v == "$" {
            n - 1
        } else if let Some(offset) = v.strip_prefix("$-") {
            (n - 1)
                .checked_sub(offset.parse()?)
                .context("index before first element")?
        } else {
            v.parse()
                .with_context(|| format!("invalid selector index '{v}'"))?
        };
        ensure!(i < n, "selector index {i} is outside 0..{}", n - 1);
        Ok(i)
    };
    let mut result = vec![];
    for item in spec.split(',') {
        let item = item.trim();
        if let Some((a, b)) = item.split_once("..") {
            let (b, step) = if let Some((b, s)) = b.split_once('(') {
                let s = s
                    .strip_suffix(')')
                    .context("selector step must end with ')'")?;
                (b, s.parse::<usize>()?)
            } else {
                (b, 1)
            };
            ensure!(step > 0, "selector step must be positive");
            let (a, b) = (index(a)?, index(b)?);
            if a <= b {
                result.extend((a..=b).step_by(step));
            } else {
                result.extend((b..=a).rev().step_by(step));
            }
        } else {
            result.push(index(item)?);
        }
    }
    ensure!(!result.is_empty(), "empty selector");
    Ok(result)
}
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum Missing {
    #[default]
    Error,
    Gap,
    Mean,
}

pub fn read(spec: &str) -> Result<Vec<Vec<f64>>> {
    read_with(spec, Missing::Error)
}

/// Select first, then handle missing values in each selected output column.
/// Inline AFNI data list samples across a line, with `|` starting a new series.
pub fn read_with(spec: &str, missing: Missing) -> Result<Vec<Vec<f64>>> {
    let transpose = spec.ends_with('\'');
    let mut name = spec.strip_suffix('\'').unwrap_or(spec);
    let mut rows = None;
    let mut cols = None;
    loop {
        let (open, is_col) = if name.ends_with(']') {
            ('[', true)
        } else if name.ends_with('}') {
            ('{', false)
        } else {
            break;
        };
        let at = name
            .rfind(open)
            .with_context(|| format!("missing '{open}' in selector"))?;
        let selection = &name[at + 1..name.len() - 1];
        if is_col {
            ensure!(cols.is_none(), "duplicate column selector");
            cols = Some(selection);
        } else {
            ensure!(rows.is_none(), "duplicate row selector");
            rows = Some(selection);
        }
        name = &name[..at];
    }
    let mut header: Vec<String> = vec![];
    let mut table: Vec<Vec<f64>> = vec![];
    if let Some(inline) = name.strip_prefix("1D:") {
        let mut series = vec![];
        for line in inline.split('|') {
            let mut values = vec![];
            for token in line
                .split(|c: char| c.is_whitespace() || c == ',')
                .filter(|t| !t.is_empty())
            {
                let (count, value) = if let Some((n, v)) = token.split_once('@') {
                    (
                        n.parse::<usize>()
                            .context("inline repeat count must be an integer")?,
                        v,
                    )
                } else {
                    (1, token)
                };
                ensure!(
                    count > 0 && count <= 1_000_000 && values.len() + count <= 1_000_000,
                    "inline data requires positive counts and at most 1000000 samples per series"
                );
                values.extend(std::iter::repeat_n(parse_value(value)?, count));
            }
            ensure!(!values.is_empty(), "empty inline series");
            if let Some(first) = series.first() {
                let first: &Vec<f64> = first;
                ensure!(first.len() == values.len(), "ragged inline series");
            }
            series.push(values);
        }
        for i in 0..series[0].len() {
            table.push(series.iter().map(|v| v[i]).collect());
        }
    } else {
        let mut text = String::new();
        if name == "-" {
            std::io::stdin().read_to_string(&mut text)?;
        } else {
            text = std::fs::read_to_string(name).with_context(|| format!("reading {name}"))?;
        }
        let tsv = name.to_ascii_lowercase().ends_with(".tsv");
        for (line_no, line) in text.lines().enumerate() {
            let line = line.split('#').next().unwrap_or("");
            if line.trim().is_empty() {
                continue;
            }
            let words: Vec<&str> = if tsv {
                line.trim_end_matches('\r')
                    .split('\t')
                    .map(str::trim)
                    .collect()
            } else {
                line.split_whitespace().collect()
            };
            if tsv && header.is_empty() && table.is_empty() {
                header = words.iter().map(|w| w.to_string()).collect();
                let mut unique = header.clone();
                unique.sort();
                unique.dedup();
                ensure!(
                    unique.len() == header.len() && header.iter().all(|s| !s.is_empty()),
                    "{name}: TSV column names must be nonempty and unique"
                );
                continue;
            }
            if let Some(first) = table.first() {
                ensure!(
                    words.len() == first.len(),
                    "{name}:{}: ragged row",
                    line_no + 1
                );
            } else if tsv {
                ensure!(words.len() == header.len(), "{name}: ragged TSV row");
            }
            // Non-numeric TSV columns remain unavailable, preserving original column indices.
            table.push(
                words
                    .iter()
                    .enumerate()
                    .map(|(c, w)| {
                        parse_value(w).or_else(|e| {
                            if tsv {
                                Ok(f64::INFINITY)
                            } else {
                                Err(e.context(format!("{name}:{} column {}", line_no + 1, c + 1)))
                            }
                        })
                    })
                    .collect::<Result<_>>()?,
            );
        }
    }
    ensure!(!table.is_empty(), "{name}: empty numeric input");
    let width = table[0].len();
    let numeric = |c: usize| table.iter().all(|r| !r[c].is_infinite());
    let ci = if let Some(c) = cols {
        let mut out = vec![];
        for part in c.split(',') {
            if let Some(i) = header.iter().position(|h| h == part.trim()) {
                out.push(i);
            } else {
                out.extend(indices(part, width)?);
            }
        }
        out
    } else {
        (0..width).filter(|&c| numeric(c)).collect()
    };
    ensure!(!ci.is_empty(), "{name}: no numeric columns");
    ensure!(
        ci.iter().all(|&c| numeric(c)),
        "{name}: selected column contains text or infinity"
    );
    let ri = if let Some(r) = rows {
        indices(r, table.len())?
    } else {
        (0..table.len()).collect()
    };
    let mut out: Vec<Vec<f64>> = ri
        .into_iter()
        .map(|r| ci.iter().map(|&c| table[r][c]).collect())
        .collect();
    if transpose {
        out = (0..out[0].len())
            .map(|c| out.iter().map(|r| r[c]).collect())
            .collect();
    }
    for c in 0..out[0].len() {
        let count = out.iter().filter(|r| r[c].is_nan()).count();
        if count == 0 {
            continue;
        }
        ensure!(
            missing != Missing::Error,
            "{name}: selected column {} has {count} missing values; use -missing gap or -missing mean",
            c + 1
        );
        if missing == Missing::Mean {
            let n = out.len() - count;
            ensure!(n > 0, "{name}: cannot impute an entirely missing column");
            let mean = out
                .iter()
                .filter(|r| r[c].is_finite())
                .map(|r| r[c] / n as f64)
                .sum::<f64>();
            ensure!(mean.is_finite(), "{name}: column mean overflow");
            for row in &mut out {
                if row[c].is_nan() {
                    row[c] = mean;
                }
            }
            eprintln!(
                "++ {name}: imputed {count} missing values in selected column {} with mean {mean}",
                c + 1
            );
        }
    }
    Ok(out)
}
fn parse_value(word: &str) -> Result<f64> {
    if ["", "na", "n/a", "nan"].contains(&word.to_ascii_lowercase().as_str()) {
        return Ok(f64::NAN);
    }
    let v: f64 = word
        .parse()
        .with_context(|| format!("expected numeric data, got '{word}'"))?;
    ensure!(!v.is_infinite(), "infinite values are unsupported");
    Ok(v)
}
pub fn column(spec: &str) -> Result<Vec<f64>> {
    let rows = read(spec)?;
    if rows[0].len() != 1 {
        bail!("{spec}: expected one column; use a '[column]' selector");
    }
    Ok(rows.into_iter().map(|r| r[0]).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn ordered_selectors() {
        assert_eq!(
            indices("$,0..4(2),3..1", 6).unwrap(),
            vec![5, 0, 2, 4, 3, 2, 1]
        );
        assert!(indices("6", 6).is_err());
        assert!(indices("0..3(0)", 6).is_err());
    }
    #[test]
    fn read_and_validate() {
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("m.1D");
        std::fs::write(&p, "# header\n1 2 3\n4 5 6 # note\n7 8 9\n").unwrap();
        assert_eq!(
            read(&format!("{}[2,0]{{2,0}}", p.display())).unwrap(),
            vec![vec![9., 7.], vec![3., 1.]]
        );
        std::fs::write(&p, "1 2\n3 NaN\n").unwrap();
        assert!(read(p.to_str().unwrap()).is_err());
        std::fs::write(&p, "1 2\n3\n").unwrap();
        assert!(read(p.to_str().unwrap()).is_err());
    }
    #[test]
    fn inline_orientation_repeats_selectors_and_transpose() {
        assert_eq!(
            read("1D: 2@0,1,2 | 4 5 6 7").unwrap(),
            vec![vec![0., 4.], vec![0., 5.], vec![1., 6.], vec![2., 7.]]
        );
        assert_eq!(
            read("1D: 1 2 3 | 4 5 6[1,0]{2,0}'").unwrap(),
            vec![vec![6., 4.], vec![3., 1.]]
        );
        assert!(read("1D: 2@1 | 1").is_err());
        assert!(read("1D: 999999999@0").is_err());
    }
    #[test]
    fn named_tsv_selection_and_explicit_missing_policies() {
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("data.tsv");
        std::fs::write(&p, "name\teffect\terror\nA\t1\t0.1\nB\tn/a\t\nC\t3\t0.3\n").unwrap();
        let name = format!("{}[effect,error]", p.display());
        assert!(read(&name).is_err());
        let m = read_with(&name, Missing::Gap).unwrap();
        assert_eq!(m.len(), 3);
        assert!(m[1][0].is_nan());
        assert!(m[1][1].is_nan());
        let m = read_with(&name, Missing::Mean).unwrap();
        assert_eq!(m[1], vec![2., 0.2]);
        assert!(read_with(&format!("{}[name]", p.display()), Missing::Gap).is_err());
        assert!(read_with("1D: 3@NaN", Missing::Mean).is_err());
        // Means use selected rows, never observations outside the plotted subset.
        assert_eq!(
            read_with("1D: 1 NaN 100{0,1}", Missing::Mean).unwrap(),
            vec![vec![1.], vec![1.]]
        );
    }
}
