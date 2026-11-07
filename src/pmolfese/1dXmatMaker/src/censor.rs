// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Combine row exclusions without changing original sample coordinates.
//! All rules are intersections of keep masks. They are applied only AFTER
//! baseline demeaning and event construction; censoring never shifts an onset.
use crate::{
    error::{Error, Result},
    input,
    model::Acquisition,
};

/// Apply per-run first/last limits and AFNI CENSORTR tokens to an owned mask.
/// Runs in text are one-based; row indices are zero-based. Out-of-range rules
/// fail instead of silently clipping. A comma starts a new independent token.
pub fn combine(
    a: &Acquisition,
    mut keep: Vec<bool>,
    first: usize,
    last: Option<usize>,
    rules: &[String],
) -> Result<Vec<bool>> {
    if keep.len() != a.nrows() {
        return Err(Error::input("censor mask length mismatch"));
    }
    for (&start, &length) in a.run_starts().iter().zip(a.run_lengths()) {
        let end = last.unwrap_or(length - 1);
        if first > end || end >= length {
            return Err(Error::input("nfirst/nlast outside a run"));
        }
        for row in 0..length {
            if row < first || row > end {
                keep[start.0 + row] = false;
            }
        }
    }
    for rule in rules {
        for token in rule
            .split(|c: char| c.is_whitespace() || c == ',')
            .filter(|s| !s.is_empty())
        {
            let (runs, selection) = if let Some((run, sel)) = token.split_once(':') {
                let runs = if run == "*" {
                    (0..a.run_lengths().len()).collect::<Vec<_>>()
                } else {
                    vec![
                        run.parse::<usize>()
                            .ok()
                            .and_then(|n| n.checked_sub(1))
                            .filter(|&n| n < a.run_lengths().len())
                            .ok_or_else(|| Error::input("invalid CENSORTR run"))?,
                    ]
                };
                (Some(runs), sel)
            } else {
                (None, token)
            };
            let selection = selection.replace('-', "..");
            if let Some(runs) = runs {
                for run in runs {
                    for i in input::indices(&selection, a.run_lengths()[run])? {
                        keep[a.run_starts()[run].0 + i] = false;
                    }
                }
            } else {
                for i in input::indices(&selection, a.nrows())? {
                    keep[i] = false;
                }
            }
        }
    }
    Ok(keep)
}
