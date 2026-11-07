// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! General linear tests map coefficients onto the FINAL design column order.
//!
//! A contrast is a small row-major matrix: each row is a linear combination
//! of fitted coefficients. Resolving names after assembly avoids confusing
//! original stimulus IDs with physical columns when baselines or AM expand.
//! This module uses Maker's own types, not afni-core's fitting machinery.
use crate::{
    design::Design,
    error::{Error, Result},
    input,
    model::normalize_label,
};
use std::path::Path;

/// Owned coefficient rows and a unique human-readable test label.
#[derive(Debug)]
pub struct Contrast {
    pub(crate) label: String,
    pub(crate) rows: Vec<Vec<f64>>,
}

impl Design {
    /// Resolve a symbolic test against this design. Ordinary `[ranges]` fill a
    /// single row; `[[ranges]]` produce one row per selected coefficient. Rust
    /// borrows the name map during parsing, then owns the finished coefficients.
    pub fn add_symbolic_contrast(&mut self, label: &str, text: &str) -> Result<()> {
        let mut rows = Vec::new();
        for line in text
            .trim()
            .strip_prefix("SYM:")
            .unwrap_or(text)
            .replace(['|', '\\'], "\n")
            .lines()
        {
            let mut ordinary = vec![0.0; self.columns().len()];
            let mut expanded: Vec<(Vec<usize>, f64)> = Vec::new();
            let mut sign = 1.0;
            let mut pending_sign = false;
            let mut found = false;
            for token in line
                .split_whitespace()
                .take_while(|t| !t.starts_with('#') && !t.starts_with("//"))
            {
                if token == "+" || token == "-" {
                    if pending_sign {
                        return Err(Error::input("consecutive GLT signs"));
                    }
                    pending_sign = true;
                    sign = if token == "-" { -1.0 } else { 1.0 };
                    continue;
                }
                found = true;
                pending_sign = false;
                let (weight, name) =
                    if let Some((w, n)) = token.split_once('*').filter(|(w, _)| !w.contains('[')) {
                        (
                            w.parse::<f64>()
                                .map_err(|_| Error::input("invalid GLT weight"))?
                                * sign,
                            n,
                        )
                    } else if let Some(n) = token.strip_prefix('-') {
                        (-sign, n)
                    } else {
                        (sign, token.strip_prefix('+').unwrap_or(token))
                    };
                let weight = (weight as f32) as f64;
                sign = 1.0;
                if !weight.is_finite() {
                    return Err(Error::input("nonfinite GLT weight"));
                }
                let (name, selector, double) = if let Some(i) = name.find('[') {
                    let double = name[i..].starts_with("[[");
                    let width = if double { 2 } else { 1 };
                    if !name.ends_with(if double { "]]" } else { "]" }) {
                        return Err(Error::input("unclosed GLT selector"));
                    }
                    (
                        &name[..i],
                        Some(&name[i + width..name.len() - width]),
                        double,
                    )
                } else {
                    (name, None, false)
                };
                if self.ranges.iter().filter(|r| r.0 == name).count() > 1 {
                    return Err(Error::input(
                        "ambiguous GLT label (Ort is reserved for polynomials)",
                    ));
                }
                let (_, first, min, count) = self
                    .ranges
                    .iter()
                    .find(|r| r.0 == name)
                    .ok_or_else(|| Error::input(format!("unknown GLT label '{name}'")))?;
                let selected = if let Some(s) = selector.filter(|s| *s != "*") {
                    input::indices(s, min + count)?
                        .into_iter()
                        .map(|i| {
                            if i < *min {
                                Err(Error::input("GLT subscript precedes minimum lag"))
                            } else {
                                Ok(first + i - min)
                            }
                        })
                        .collect::<Result<Vec<_>>>()?
                } else {
                    (*first..first + count).collect()
                };
                if double {
                    expanded.push((selected, weight));
                } else {
                    for c in selected {
                        ordinary[c] = weight;
                    }
                }
            }
            if pending_sign {
                return Err(Error::input("trailing GLT sign"));
            }
            if !found {
                continue;
            }
            if expanded.is_empty() {
                rows.push(ordinary);
            } else {
                let n = expanded[0].0.len();
                if expanded.iter().any(|(s, _)| s.len() != n) {
                    return Err(Error::input(
                        "double-bracket GLT ranges must have equal lengths",
                    ));
                }
                for i in 0..n {
                    // Expanded coefficients are assigned first. AFNI then
                    // fills every remaining zero from the ordinary row; an
                    // explicit expanded zero therefore does not erase it.
                    let mut row = vec![0.0; ordinary.len()];
                    for (selected, w) in &expanded {
                        row[selected[i]] = *w;
                    }
                    for (value, &fallback) in row.iter_mut().zip(&ordinary) {
                        if *value == 0.0 {
                            *value = fallback;
                        }
                    }
                    rows.push(row);
                }
            }
        }
        self.push_contrast(label, rows)
    }
    /// Read numeric tests as exactly nrows by all design columns, including
    /// baseline zeros. No implicit baseline padding or truncated rows is used.
    pub fn add_numeric_contrast(&mut self, label: &str, path: &Path, nrows: usize) -> Result<()> {
        let text = input::rows(path, false)?;
        if text.len() != nrows {
            return Err(Error::input("numeric GLT row count mismatch"));
        }
        let rows = text
            .iter()
            .map(|r| {
                r.iter()
                    .map(|v| {
                        v.parse::<f64>()
                            .map_err(|_| Error::input("invalid numeric GLT"))
                    })
                    .collect()
            })
            .collect::<Result<Vec<Vec<f64>>>>()?;
        self.push_contrast(label, rows)
    }
    /// Central validation keeps symbolic and numeric tests under one contract.
    fn push_contrast(&mut self, label: &str, mut rows: Vec<Vec<f64>>) -> Result<()> {
        // The reference GLT reader stores float coefficients before widening
        // into its matrix. Round once here, including for numeric input files.
        for value in rows.iter_mut().flatten() {
            *value = (*value as f32) as f64;
        }
        let label = normalize_label(label)?;
        if rows.is_empty()
            || rows
                .iter()
                .any(|r| r.len() != self.columns().len() || r.iter().any(|v| !v.is_finite()))
        {
            return Err(Error::input(
                "GLT must contain finite, full-width coefficient rows",
            ));
        }
        if self.contrasts.iter().any(|g| g.label == label) {
            return Err(Error::input("duplicate GLT label"));
        }
        self.contrasts.push(Contrast { label, rows });
        Ok(())
    }
}
