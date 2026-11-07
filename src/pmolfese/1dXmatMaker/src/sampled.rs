// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Sampled responses and the pinned, real GLMsingle 20-HRF library.
//!
//! Read Curve first (uniform samples and linear evaluation), then glmsingle
//! (box summation, PCHIP resampling, normalization). No AFNI or Python call is
//! made here: the embedded TSV is the upstream numeric library, not twenty
//! invented analytic approximations. Independent tests run upstream gethrf.py.
//!
//! GLMsingle data/algorithm provenance: cvnlab/GLMsingle commit
//! 082851f9db87ff30c49370f27365ac1112f18ff7, Kendrick Kay, BSD-3-Clause.
//! The complete license is retained in data/glmsingle/LICENSE.
use crate::error::{Error, Result};

/// Pinned library revision; changes require an explicit new fixture capture.
pub const GLMSINGLE_REVISION: &str = "082851f9db87ff30c49370f27365ac1112f18ff7";
/// Owned, uniformly sampled, already prepared response. Duration has already
/// been handled: evaluating this curve must NEVER convolve a second time.
#[derive(Debug, Clone)]
pub struct Curve {
    values: Vec<f64>,
    dt: f64,
    duration: f64,
}
impl Curve {
    /// Validate finite samples, positive spacing and nonnegative duration.
    /// The first sample is at onset (zero seconds). Values outside the sampled
    /// support are zero; within it, fractional onsets use linear interpolation.
    pub fn new(values: Vec<f64>, dt: f64, duration: f64) -> Result<Self> {
        if values.len() < 2
            || values.len() > 1_000_000
            || values.iter().any(|v| !v.is_finite())
            || !dt.is_finite()
            || dt <= 0.0
            || !duration.is_finite()
            || duration < 0.0
            || !((values.len() - 1) as f64 * dt).is_finite()
        {
            return Err(Error::input(
                "sampled HRF needs 2..1000000 finite values, positive dt and finite duration",
            ));
        }
        Ok(Self {
            values,
            dt,
            duration,
        })
    }
    /// Borrow waveform samples. No copy or external allocation is returned.
    pub fn values(&self) -> &[f64] {
        &self.values
    }
    /// Grid spacing in seconds.
    pub fn dt(&self) -> f64 {
        self.dt
    }
    /// Fixed event duration already represented by this curve.
    pub fn duration(&self) -> f64 {
        self.duration
    }
    /// Last sampled time, relative to onset, in seconds.
    pub fn end(&self) -> f64 {
        (self.values.len() - 1) as f64 * self.dt
    }
    /// Linear sampling; exact knots preserve original values and tails are zero.
    pub fn sample(&self, time: f64) -> Option<f64> {
        if !time.is_finite() {
            return None;
        }
        if time < 0.0 || time > self.end() {
            return Some(0.0);
        }
        let x = time / self.dt;
        let i = (x.floor() as usize).min(self.values.len() - 1);
        if i + 1 == self.values.len() {
            return Some(self.values[i]);
        }
        Some(self.values[i] + (x - i as f64) * (self.values[i + 1] - self.values[i]))
    }
    /// Normalize the positive peak to one. Negative-only or zero responses
    /// cannot use this convention; callers may preserve their original units.
    pub fn peak_normalized(mut self) -> Result<Self> {
        let peak = self
            .values
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        if peak <= 0.0 {
            return Err(Error::input(
                "peak normalization requires a positive HRF sample",
            ));
        }
        for value in &mut self.values {
            *value /= peak;
        }
        // Very small positive peaks can overflow negative samples on division.
        // Reuse the constructor so normalization never exports nonfinite data.
        Self::new(self.values, self.dt, self.duration)
    }
    /// Integrate an impulse response against a unit-height box of duration d.
    /// The integral of a piecewise-linear segment is a trapezoid, evaluated
    /// exactly, including fractional duration endpoints. d=0 means impulse,
    /// not a zero-area box. Unlike GLMsingle's 0.1-s response library, custom
    /// impulse samples are integrated in seconds and retain amplitude units.
    pub fn convolve_box(&self, d: f64) -> Result<Self> {
        if !d.is_finite() || !(0.0..=999.0).contains(&d) {
            return Err(Error::input("custom duration must be in 0..999 seconds"));
        }
        if d == 0.0 {
            return Self::new(self.values.clone(), self.dt, 0.0);
        }
        // Check in floating point BEFORE converting or adding one. Rust casts
        // enormous floats to usize::MAX; adding one would overflow that value.
        let intervals = ((self.end() + d) / self.dt).ceil();
        if !intervals.is_finite() || intervals > 999_999.0 {
            return Err(Error::input("convolved HRF exceeds sample limit"));
        }
        let count = intervals as usize + 1;
        let mut area = vec![0.0; self.values.len()];
        for i in 1..area.len() {
            area[i] = area[i - 1] + self.dt * (self.values[i - 1] + self.values[i]) / 2.0;
        }
        let integral = |time: f64| {
            if time <= 0.0 {
                return 0.0;
            }
            if time >= self.end() {
                return area[area.len() - 1];
            }
            let i = (time / self.dt).floor() as usize;
            let t = time - i as f64 * self.dt;
            area[i]
                + self.values[i] * t
                + (self.values[i + 1] - self.values[i]) * t * t / (2.0 * self.dt)
        };
        Self::new(
            (0..count)
                .map(|i| {
                    let t = i as f64 * self.dt;
                    integral(t) - integral(t - d)
                })
                .collect(),
            self.dt,
            d,
        )
    }
}

/// Prepare all twenty upstream HRFs for a fixed duration and acquisition TR.
/// Upstream getcanonicalhrflibrary first globally peak-normalizes then returns
/// float32. Its fit() additionally normalizes each candidate's peak; `per_peak`
/// selects that second step. Keep these distinct for faithful oracle tests.
///
/// The upstream samples represent responses to a 0.1-second stimulus. Duration
/// is rounded half away from zero into that many 0.1-s boxes, at least one.
/// PCHIP includes upstream's final-grid extrapolation before float32 conversion.
pub fn glmsingle(duration: f64, tr: f64, per_peak: bool) -> Result<Vec<Curve>> {
    if !duration.is_finite() || !(0.0..=999.0).contains(&duration) || !tr.is_finite() || tr <= 0.0 {
        return Err(Error::input(
            "GLMsingle requires duration 0..999 and positive finite TR",
        ));
    }
    let text = include_str!("../data/glmsingle/getcanonicalhrflibrary.tsv");
    let mut columns = vec![Vec::new(); 20];
    for line in text.lines() {
        let row = line
            .split_whitespace()
            .map(|v| {
                v.parse::<f64>()
                    .map_err(|_| Error::input("invalid embedded HRF number"))
            })
            .collect::<Result<Vec<_>>>()?;
        if row.len() != 20 {
            return Err(Error::input(
                "embedded GLMsingle library must have twenty columns",
            ));
        }
        for (column, value) in columns.iter_mut().zip(row) {
            column.push(value);
        }
    }
    if columns[0].len() != 501 {
        return Err(Error::input(
            "embedded GLMsingle library must have 501 rows",
        ));
    }
    let width = (duration / 0.1).round().max(1.0) as usize;
    let length = 501 + width - 1;
    let stop = (length as f64 * 0.1).ceil();
    let count = (stop / tr).ceil() as usize;
    if !(2..=1_000_000).contains(&count) {
        return Err(Error::input(
            "resampled GLMsingle grid needs 2..1000000 samples",
        ));
    }
    let mut sampled = Vec::new();
    let mut peak = f64::NEG_INFINITY;
    for original in columns {
        // Sum boxes in source order, like numpy.convolve. A sliding sum would
        // accumulate different cancellation error in the small negative tail.
        let mut convolved = vec![0.0; length];
        for (i, &v) in original.iter().enumerate() {
            for target in &mut convolved[i..i + width] {
                *target += v;
            }
        }
        let slopes = pchip_slopes(&convolved, 0.1);
        let row = (0..count)
            .map(|i| pchip_value(&convolved, &slopes, 0.1, i as f64 * tr))
            .collect::<Vec<_>>();
        peak = peak.max(row.iter().copied().fold(f64::NEG_INFINITY, f64::max));
        sampled.push(row);
    }
    if !peak.is_finite() || peak <= 0.0 {
        return Err(Error::input("invalid GLMsingle normalization peak"));
    }
    let mut result = Vec::new();
    for mut row in sampled {
        for v in &mut row {
            *v = (*v / peak) as f32 as f64;
        }
        if per_peak {
            let local = row.iter().copied().fold(f64::NEG_INFINITY, f64::max) as f32;
            for v in &mut row {
                *v = (*v as f32 / local) as f64;
            }
        }
        result.push(Curve::new(row, tr, duration)?);
    }
    Ok(result)
}

/// Fritsch–Butland monotone cubic slopes on a uniform grid. Opposite-sign
/// secants give zero slope at extrema; same-sign neighbors use a harmonic
/// mean. One-sided endpoint estimates are clipped to avoid overshoot. This is
/// the PCHIP convention used by scipy.interpolate.PchipInterpolator upstream.
fn pchip_slopes(y: &[f64], dt: f64) -> Vec<f64> {
    let d = y.windows(2).map(|w| (w[1] - w[0]) / dt).collect::<Vec<_>>();
    let n = y.len();
    let mut m = vec![0.0; n];
    for i in 1..n - 1 {
        if d[i - 1] != 0.0 && d[i] != 0.0 && d[i - 1].signum() == d[i].signum() {
            m[i] = 2.0 / (1.0 / d[i - 1] + 1.0 / d[i]);
        }
    }
    let endpoint = |a: f64, b: f64| {
        let value = (3.0 * a - b) / 2.0;
        if value.signum() != a.signum() {
            0.0
        } else if a.signum() != b.signum() && value.abs() > 3.0 * a.abs() {
            3.0 * a
        } else {
            value
        }
    };
    m[0] = endpoint(d[0], d[1]);
    m[n - 1] = endpoint(d[n - 2], d[n - 3]);
    m
}
/// Cubic Hermite evaluation, including extrapolation of the last interval.
/// Only library PREPARATION extrapolates; the final Curve has zero exterior.
fn pchip_value(y: &[f64], m: &[f64], dt: f64, t: f64) -> f64 {
    let i = ((t / dt).floor() as usize).min(y.len() - 2);
    let u = (t - i as f64 * dt) / dt;
    (2.0 * u * u * u - 3.0 * u * u + 1.0) * y[i]
        + (u * u * u - 2.0 * u * u + u) * dt * m[i]
        + (-2.0 * u * u * u + 3.0 * u * u) * y[i + 1]
        + (u * u * u - u * u) * dt * m[i + 1]
}
