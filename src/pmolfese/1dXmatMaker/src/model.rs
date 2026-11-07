// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Validated acquisition and row identities, independent of file parsing.
//!
//! Run lengths are original acquisition samples, TR is seconds, and every
//! index here is zero-based. Small wrapper types prevent passing a retained
//! row or a run number where a full-timeline row is expected. Private fields
//! keep constructors' invariants intact for every library caller.
//!
//! Rust reading note: `struct FullRow(pub usize)` is a distinct type, not a
//! typedef. FullRow(3) and RetainedRow(3) cannot be exchanged accidentally.
//! `.0` extracts their integer when an index is finally needed. `Copy` is safe
//! for these tiny integers; Acquisition deliberately owns its vectors instead.
//! No afni-core call occurs here: these are Maker's local acquisition types.

use crate::error::{Error, Result, reserve};

/// Zero-based run identity, unrelated to full-timeline row offsets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RunIndex(pub usize);
/// Zero-based index on the original, uncensored acquisition timeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FullRow(pub usize);
/// Zero-based position among retained rows, after censoring.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RetainedRow(pub usize);
/// Zero-based matrix column; this is not a fitted beta-volume index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ColumnIndex(pub usize);

/// Explicit acquisition settings, validated once before matrix allocation.
#[derive(Debug, Clone)]
pub struct Acquisition {
    lengths: Vec<usize>,
    starts: Vec<FullRow>,
    nrows: usize,
    tr: f64,
    polort: i32,
    demean: bool,
    power: bool,
}

impl Acquisition {
    /// Own run lengths and validate TR, degree, and cumulative row counts.
    ///
    /// `polort=-1` omits polynomials AND disables demeaning, even if explicitly
    /// requested (3dDeconvolve.c line 3897 at the M0 reference). Otherwise
    /// `demean=None` enables demeaning. Empty/zero-length runs, nonfinite TR, one-sample
    /// nonconstant trends, and dimensions exceeding AFNI's i32 indices fail.
    pub fn new(lengths: Vec<usize>, tr: f64, polort: i32, demean: Option<bool>) -> Result<Self> {
        if lengths.is_empty() || !tr.is_finite() || tr <= 0.0 {
            return Err(Error::input(
                "require nonempty run lengths and finite positive TR",
            ));
        }
        if polort < -1 {
            return Err(Error::input(
                "library polort must be -1 or a nonnegative integer; the CLI resolves A",
            ));
        }
        let mut starts = Vec::new();
        reserve(&mut starts, lengths.len())?;
        let mut nrows = 0usize;
        for (run, &length) in lengths.iter().enumerate() {
            if length == 0 || (length == 1 && polort > 0) {
                return Err(Error::input(format!(
                    "run {}: length must be positive; nonconstant polynomials need at least two samples",
                    run + 1
                )));
            }
            starts.push(FullRow(nrows));
            // checked_add returns Option: None means integer overflow. The
            // filter also enforces AFNI's index width before any allocation.
            nrows = nrows
                .checked_add(length)
                .filter(|&n| n <= i32::MAX as usize)
                .ok_or_else(|| Error::input("full row count exceeds AFNI's 32-bit index range"))?;
            if !(length as f64 * tr).is_finite() {
                return Err(Error::input(format!(
                    "run {}: duration overflows finite seconds",
                    run + 1
                )));
            }
        }
        Ok(Self {
            lengths,
            starts,
            nrows,
            tr,
            polort,
            demean: polort >= 0 && demean.unwrap_or(true),
            power: false,
        })
    }

    /// Select legacy run-local powers t^k instead of Legendre Pk. This consumes
    /// the acquisition and returns it, so no other owner can observe a partial
    /// change. Power coordinates are sample indices, not normalized seconds.
    pub fn with_power_baseline(mut self, power: bool) -> Self {
        self.power = power;
        self
    }
    /// Whether baseline columns use legacy powers of run-local sample index.
    pub fn power_baseline(&self) -> bool {
        self.power
    }

    /// Borrow original run lengths. The slice is valid while `self` is alive.
    pub fn run_lengths(&self) -> &[usize] {
        &self.lengths
    }
    /// Borrow starts on the full timeline; censoring never changes them.
    pub fn run_starts(&self) -> &[FullRow] {
        &self.starts
    }
    /// Number of full-timeline samples, before any row removal.
    pub fn nrows(&self) -> usize {
        self.nrows
    }
    /// Common acquisition interval in seconds.
    pub fn tr(&self) -> f64 {
        self.tr
    }
    /// Polynomial degree, or -1 when there are no polynomial columns.
    pub fn polort(&self) -> i32 {
        self.polort
    }
    /// Whether to remove means from nonconstant trends and nuisance vectors.
    pub fn demean(&self) -> bool {
        self.demean
    }
    /// Return one run's interval, or None for an out-of-range run identity.
    pub fn run_range(&self, run: RunIndex) -> Option<std::ops::Range<usize>> {
        let start = self.starts.get(run.0)?.0;
        Some(start..start + self.lengths[run.0])
    }
}

/// Normalize spaces as AFNI does, while refusing ambiguous metadata labels.
///
/// AFNI's 3dDeconvolve column formatter retains 60 bytes before appending #0
/// (source column_metadata is 64 bytes). Require the *expanded* label to fit
/// rather than silently truncating it; callers check ortvec's [column] suffix.
pub(crate) fn normalize_label(label: &str) -> Result<String> {
    if label.is_empty() || label.chars().any(|c| c.is_control() || c == ';') {
        return Err(Error::input(
            "labels must be nonempty and contain no semicolons or control characters",
        ));
    }
    let normalized: String = label
        .chars()
        .map(|c| if c.is_whitespace() { '_' } else { c })
        .collect();
    if normalized.len() > 60 {
        return Err(Error::input(format!(
            "label '{label}' exceeds the 60-byte AFNI compatibility limit"
        )));
    }
    Ok(normalized)
}
