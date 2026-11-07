// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Plain numeric tables: finite f64 values, whitespace, and # comments only.
//!
//! A table owns one column-major buffer. The text reader first collects rows,
//! then transposes once so design assembly can borrow complete columns.
//! input.rs applies selectors and inline syntax first. Text-to-zero repair is
//! never allowed. Errors identify the source and selected/transformed line.
//!
//! File I/O boundary: input.rs uses the standard library, not afni-io. The
//! Cursor below adapts selected text to the same buffered-parser interface used
//! by tests. read_table returns an owned NumericTable, so no
//! numerical column points into a temporary text buffer after reading finishes.
//! If afni-io is introduced later, document its exact ownership, dimensions,
//! precision and error conversions at the call site rather than hiding them.

use crate::error::{Error, Result, reserve};
use std::{io::BufRead, path::Path};

/// A finite rectangular numeric table with private validated storage.
#[derive(Debug, Clone)]
pub struct NumericTable {
    nrows: usize,
    ncols: usize,
    values: Vec<f64>,
}

impl NumericTable {
    /// Own a column-major matrix after checking its shape and finite values.
    ///
    /// `values[c*nrows+r]` must be row r of column c. Zero dimensions and
    /// mismatched/overflowing dimensions fail before values can be indexed.
    pub fn new(nrows: usize, ncols: usize, values: Vec<f64>) -> Result<Self> {
        if nrows == 0 || ncols == 0 || nrows.checked_mul(ncols) != Some(values.len()) {
            return Err(Error::input(
                "numeric table must be nonempty and rectangular",
            ));
        }
        if values.iter().any(|v| !v.is_finite()) {
            return Err(Error::input("numeric table contains a nonfinite value"));
        }
        Ok(Self {
            nrows,
            ncols,
            values,
        })
    }

    /// Number of time samples in this full-timeline table.
    pub fn nrows(&self) -> usize {
        self.nrows
    }
    /// Number of supplied vectors.
    pub fn ncols(&self) -> usize {
        self.ncols
    }
    /// Borrow one contiguous column; None means the column does not exist.
    pub fn column(&self, column: usize) -> Option<&[f64]> {
        if column >= self.ncols {
            return None;
        }
        let start = column * self.nrows;
        Some(&self.values[start..start + self.nrows])
    }
}

/// Read a numeric file requiring exactly `expected_rows` original samples.
///
/// Length mismatches, ragged rows, nonfinite tokens, and unsupported selector
/// syntax are input errors. Opening/reading a file is an I/O error. AFNI text
/// selectors are applied by input.rs; shell variables, path globs and ~ are
/// never expanded by this reader.
pub fn read_table(path: &Path, expected_rows: usize) -> Result<NumericTable> {
    let rows = crate::input::rows(path, false)?;
    let text = rows
        .iter()
        .map(|r| r.join(" "))
        .collect::<Vec<_>>()
        .join("\n");
    read_numeric(
        std::io::Cursor::new(text),
        &path.to_string_lossy(),
        expected_rows,
    )
}

/// Parse a buffered stream; isolated from file opening for in-memory tests.
fn read_numeric(reader: impl BufRead, source: &str, expected_rows: usize) -> Result<NumericTable> {
    let mut rows = Vec::new();
    let mut nrows = 0usize;
    let mut ncols = None;
    for (line_number, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| {
            Error::resource(format!(
                "{source}: cannot read line {}: {e}",
                line_number + 1
            ))
        })?;
        let content = line.split('#').next().unwrap_or("").trim();
        if content.is_empty() {
            continue;
        }
        if nrows >= expected_rows {
            return Err(Error::input(format!(
                "{source}: line {} exceeds expected {expected_rows} full rows",
                line_number + 1
            )));
        }
        let tokens: Vec<_> = content.split_whitespace().collect();
        let width = *ncols.get_or_insert(tokens.len());
        if tokens.len() != width {
            return Err(Error::input(format!(
                "{source}: line {} has {} columns, expected {width}",
                line_number + 1,
                tokens.len()
            )));
        }
        reserve(&mut rows, width)?;
        for (column, token) in tokens.iter().enumerate() {
            let value: f64 = token.parse().map_err(|_| {
                Error::input(format!(
                    "{source}: line {}, column {}: invalid number '{token}'",
                    line_number + 1,
                    column + 1
                ))
            })?;
            if !value.is_finite() {
                return Err(Error::input(format!(
                    "{source}: line {}, column {}: nonfinite value",
                    line_number + 1,
                    column + 1
                )));
            }
            rows.push(value);
        }
        nrows += 1;
    }
    if nrows != expected_rows || nrows == 0 {
        return Err(Error::input(format!(
            "{source}: has {nrows} numeric rows, expected {expected_rows}"
        )));
    }
    let width = ncols.unwrap_or(0);
    let mut values = Vec::new();
    reserve(&mut values, rows.len())?;
    // Input arrived row by row (r*width+c). Store each whole column together
    // instead, so later demeaning and regression construction use one slice.
    // This is a storage transpose, not a change in the scientific matrix.
    for column in 0..width {
        for row in 0..nrows {
            values.push(rows[row * width + column]);
        }
    }
    NumericTable::new(nrows, width, values)
}

/// Read a single 0/1 censor column: true retains the original row.
pub fn read_censor(path: &Path, expected_rows: usize) -> Result<Vec<bool>> {
    let table = read_table(path, expected_rows)?;
    if table.ncols != 1 {
        return Err(Error::input(format!(
            "{}: censor file must have one column",
            path.display()
        )));
    }
    let mut keep = Vec::new();
    reserve(&mut keep, expected_rows)?;
    for (row, &value) in table.values.iter().enumerate() {
        if value != 0.0 && value != 1.0 {
            return Err(Error::input(format!(
                "{}: censor row {} must be 0 or 1",
                path.display(),
                row + 1
            )));
        }
        keep.push(value == 1.0);
    }
    Ok(keep)
}
