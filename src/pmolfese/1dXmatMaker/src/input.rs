// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Text-only AFNI input syntax shared by tables, timing and GLT selectors.
//!
//! No afni-io call occurs here. The standard library reads UTF-8 and this
//! module interprets the documented subset of AFNI's 1D syntax. Selection is
//! done BEFORE shape validation; a selector can intentionally reorder/repeat
//! rows or columns. It cannot turn bad text into a numeric zero.
use crate::error::{Error, Result};
use std::{fs, path::Path};

/// Expand an inclusive integer list: commas, a..b, a..b(step), $, or *.
/// `length` defines the valid zero-based index range. A Vec owns the result;
/// callers can then select from any borrowed table without retaining text.
pub fn indices(text: &str, length: usize) -> Result<Vec<usize>> {
    if length == 0 {
        return Err(Error::input("selector applied to an empty axis"));
    }
    let number = |s: &str| -> Result<usize> {
        let n = if let Some(rest) = s.trim().strip_prefix('$') {
            if rest.is_empty() {
                length - 1
            } else {
                let minus = rest
                    .strip_prefix('-')
                    .ok_or_else(|| Error::input("invalid $ offset"))?;
                (length - 1)
                    .checked_sub(
                        minus
                            .parse()
                            .map_err(|_| Error::input("invalid $ offset"))?,
                    )
                    .ok_or_else(|| Error::input("$ offset exceeds axis"))?
            }
        } else {
            s.trim()
                .parse()
                .map_err(|_| Error::input("invalid selector index"))?
        };
        if n >= length {
            return Err(Error::input("selector index exceeds axis"));
        }
        Ok(n)
    };
    if length > 16_777_216 {
        return Err(Error::input("selector axis exceeds resource limit"));
    }
    let mut out = Vec::new();
    for part in text.split(',') {
        if part.trim() == "*" {
            if out.len() + length > 16_777_216 {
                return Err(Error::input("selector is too large"));
            }
            out.extend(0..length);
            continue;
        }
        let (range, step) = if let Some((r, s)) = part.split_once('(') {
            let s = s
                .strip_suffix(')')
                .ok_or_else(|| Error::input("unterminated selector step"))?;
            (
                r,
                s.parse::<usize>()
                    .map_err(|_| Error::input("invalid selector step"))?,
            )
        } else {
            (part, 1)
        };
        if step == 0 {
            return Err(Error::input("selector step must be positive"));
        }
        if let Some((a, b)) = range.split_once("..") {
            let (a, b) = (number(a)?, number(b)?);
            if a <= b {
                out.extend((a..=b).step_by(step));
            } else {
                out.extend((b..=a).rev().step_by(step));
            }
        } else {
            out.push(number(range)?);
        }
        if out.len() > 16_777_216 {
            return Err(Error::input("selector is too large"));
        }
    }
    if out.is_empty() {
        return Err(Error::input("empty selector"));
    }
    Ok(out)
}

/// Read token rows, applying optional `[columns]`, `{rows}`, and trailing transpose.
/// Inline 1D is ONE column; commas/semicolons separate values and n@value repeats.
/// Timing uses the file orientation, while inline timing is a single row.
pub fn rows(spec: &Path, timing: bool) -> Result<Vec<Vec<String>>> {
    let name = spec
        .to_str()
        .ok_or_else(|| Error::input("input path must be UTF-8"))?;
    let (name, transpose) = name.strip_suffix('\'').map_or((name, false), |n| (n, true));
    let mut base = name;
    let mut columns = None;
    let mut selected_rows = None;
    // Peel selectors from the right so filenames and the inline payload are
    // left intact. Duplicate selectors are rejected rather than reordered.
    loop {
        let (open, close) = if base.ends_with(']') {
            ('[', ']')
        } else if base.ends_with('}') {
            ('{', '}')
        } else {
            break;
        };
        let i = base
            .rfind(open)
            .ok_or_else(|| Error::input("unmatched input selector"))?;
        let content = &base[i + 1..base.len() - close.len_utf8()];
        let slot = if open == '[' {
            &mut columns
        } else {
            &mut selected_rows
        };
        if slot.replace(content).is_some() {
            return Err(Error::input("duplicate input selector"));
        }
        base = &base[..i];
    }
    let mut rows = if let Some(payload) = base.strip_prefix("1D:") {
        let mut cols = Vec::new();
        for column in [payload] {
            let mut values = Vec::new();
            for token in column
                .split(|c: char| c.is_whitespace() || c == ';' || (!timing && c == ','))
                .filter(|s| !s.is_empty())
            {
                let (count, value) = if let Some((n, v)) = token.split_once('@') {
                    (
                        n.parse::<usize>()
                            .map_err(|_| Error::input("invalid inline repetition"))?,
                        v,
                    )
                } else {
                    (1, token)
                };
                if count == 0
                    || values
                        .len()
                        .checked_add(count)
                        .is_none_or(|n| n > 16_777_216)
                {
                    return Err(Error::input("inline repetition is empty or too large"));
                }
                values.extend(std::iter::repeat_n(value.to_owned(), count));
            }
            cols.push(values);
        }
        if timing {
            if cols.len() != 1 {
                return Err(Error::input(
                    "inline timing requires a single component row",
                ));
            }
            cols
        } else {
            transpose_rows(cols)?
        }
    } else {
        if base == "-" {
            return Err(Error::input(
                "stdin is unsupported; use a file or 1D: literal",
            ));
        }
        let text =
            fs::read_to_string(base).map_err(|e| Error::io("cannot read", Path::new(base), e))?;
        text.lines()
            .filter_map(|line| {
                let line = line.split('#').next().unwrap_or("").trim();
                if line.is_empty() {
                    None
                } else {
                    Some(line.split_whitespace().map(str::to_owned).collect())
                }
            })
            .collect()
    };
    if rows.is_empty() {
        return Err(Error::input("input has no data rows"));
    }
    if let Some(selection) = columns {
        // Ragged timing is allowed, but selecting a missing token is an error.
        let chosen = indices(selection, rows[0].len())?;
        rows = rows
            .iter()
            .map(|row| {
                chosen
                    .iter()
                    .map(|&c| {
                        row.get(c)
                            .cloned()
                            .ok_or_else(|| Error::input("column selector exceeds a ragged row"))
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
    }
    if let Some(selection) = selected_rows {
        rows = indices(selection, rows.len())?
            .iter()
            .map(|&r| rows[r].clone())
            .collect();
    }
    if transpose {
        rows = transpose_rows(rows)?;
    }
    Ok(rows)
}

/// Rectangular transpose of owned token strings. This is an input transform,
/// distinct from column-major numeric storage in table.rs.
fn transpose_rows(rows: Vec<Vec<String>>) -> Result<Vec<Vec<String>>> {
    if rows.is_empty() || rows[0].is_empty() || rows.iter().any(|r| r.len() != rows[0].len()) {
        return Err(Error::input(
            "transpose/inline columns require a nonempty rectangle",
        ));
    }
    Ok((0..rows[0].len())
        .map(|c| rows.iter().map(|r| r[c].clone()).collect())
        .collect())
}
