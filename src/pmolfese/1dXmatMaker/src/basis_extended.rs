// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Additional response families and AFNI-compatible duration convolution.
//!
//! Scientific behavior is adapted from 3dDeconvolve.c, cs_gamfit.c and
//! powell_int.c at the pinned revision in lib.rs. Preserve the attribution in
//! basis.rs: AFNI / B. Douglas Ward / Medical College of Wisconsin, CC BY 4.0.
//! This is Rust code, not a call into AFNI. AFNI executables are only test oracles.
//!
//! `Extended` is an enum: each variant owns just the parameters needed by one
//! family. A `Table` owns its precomputed samples. Borrowing `&self` during
//! sampling avoids copying those arrays for every scan. Durations use AFNI's
//! historical float grid; changing this to a more accurate integral would
//! change the design matrix, particularly derivative columns.

use crate::{
    basis::block_integral,
    error::{Error, Result},
    expression::Expression,
};

const DT: f32 = 0.01;
// Explicit resource limits prevent a malformed model from requesting an
// unbounded allocation or quadrature loop. They are errors, never truncation.
const MAX_GRID: usize = 200_000;
const MAX_CONVOLUTION_TERMS: usize = 100_000_000;

#[derive(Debug)]
pub(crate) enum Extended {
    Gamma {
        p: f32,
        q: f32,
        end: f32,
    },
    TwoGamma {
        parameters: [f32; 5],
        end: f32,
        scale: f32,
    },
    Block {
        duration: f32,
        peak: f32,
        power: i32,
        unit: bool,
        dynamic: bool,
    },
    Grid {
        start: f32,
        end: f32,
        count: usize,
        family: GridFamily,
    },
    Spmg(usize),
    Mion(f32),
    Table {
        duration: f32,
        curves: Vec<Vec<f32>>,
    },
    Expr {
        start: f32,
        end: f32,
        programs: Vec<Expression>,
    },
}

#[derive(Debug, Clone, Copy)]
pub(crate) enum GridFamily {
    Spline { zero: bool },
    Poly,
    Sin,
    Trig,
}

/// Return None only for families handled by basis.rs. A recognized family
/// with bad arguments returns Some(Err): it must not fall through to another
/// model. This nested return type distinguishes “not mine” from “invalid”.
pub(crate) fn parse(formula: &str) -> Option<Result<Extended>> {
    let name = formula.split('(').next().unwrap_or("");
    let handled = matches!(
        name,
        "GAMpw"
            | "TWOGAM"
            | "TWOGAMpw"
            | "IGFUN"
            | "IGFUN4"
            | "IGFUN5"
            | "UBLOCK"
            | "UBLOCK4"
            | "UBLOCK5"
            | "dmBLOCK"
            | "dmBLOCK4"
            | "dmBLOCK5"
            | "dmUBLOCK"
            | "dmUBLOCK4"
            | "dmUBLOCK5"
            | "CSPLIN"
            | "CSPLINzero"
            | "POLY"
            | "SIN"
            | "TRIG"
            | "SPMG"
            | "SPMG3"
            | "MION"
            | "MIONN"
            | "WAV"
            | "EXPR"
    ) || (matches!(name, "SPMG1" | "SPMG2") && formula.contains('('))
        || (name == "GAM" && formula.matches(',').count() == 2);
    if !handled {
        return None;
    }
    Some(parse_known(formula, name))
}

fn invalid(formula: &str) -> Error {
    Error::input(format!(
        "invalid basis '{formula}'; see -help for parameter bounds"
    ))
}

/// Parse finite float parameters without C scanf's silent trailing-text or
/// partial-argument acceptance. WAV alone permits empty fields for defaults.
fn numbers(formula: &str, allow_empty: bool) -> Result<Vec<Option<f32>>> {
    let Some((_, rest)) = formula.split_once('(') else {
        return Ok(Vec::new());
    };
    let body = rest.strip_suffix(')').ok_or_else(|| invalid(formula))?;
    if body.is_empty() {
        return if allow_empty {
            Ok(Vec::new())
        } else {
            Err(invalid(formula))
        };
    }
    body.split(',')
        .map(|s| {
            if allow_empty && s.trim().is_empty() {
                return Ok(None);
            }
            let x: f32 = s.trim().parse().map_err(|_| invalid(formula))?;
            if !x.is_finite() {
                return Err(invalid(formula));
            }
            Ok(Some(x))
        })
        .collect()
}

/// Dispatch only recognized names. Parameter validation happens before tables
/// allocate memory; the returned enum owns everything sampling will need.
fn parse_known(formula: &str, name: &str) -> Result<Extended> {
    if name == "EXPR" {
        return parse_expression(formula);
    }
    let optional = numbers(formula, name == "WAV")?;
    let p: Vec<f32> = optional.iter().map(|v| v.unwrap_or(0.0)).collect();
    let bad = || invalid(formula);
    match name {
        "GAM" | "GAMpw" => {
            if p.is_empty() && !formula.contains('(') {
                return Ok(Extended::Gamma {
                    p: 8.6,
                    q: 0.547,
                    end: 11.1,
                });
            }
            if !(p.len() == 2 || p.len() == 3) {
                return Err(bad());
            }
            let (mut shape, mut scale) = (p[0], p[1]);
            let duration = p.get(2).copied();
            if duration.is_some_and(|d| d < 0.0) {
                return Err(bad());
            }
            if name == "GAMpw" {
                (shape, scale) = peak_width(shape, scale)?;
            } else if duration.is_some() {
                // In AFNI, explicit duration lets nonpositive p/q request the
                // defaults. Without duration those values are invalid.
                if shape <= 0.0 {
                    shape = 8.6;
                }
                if scale <= 0.0 {
                    scale = 0.547;
                }
            }
            if shape <= 0.0 || scale <= 0.0 {
                return Err(bad());
            }
            let end = gamma_end(shape, scale, if duration.is_some() { 5.0 } else { 9.9 })?;
            if let Some(d) = duration {
                // GAM(…,0) intentionally uses TWO microtime lags, whereas
                // SPMG(0) uses one. This surprising difference is in AFNI.
                let d = if d == 0.0 { DT } else { d };
                table(d, 1, end, |_, t| gamma(t, shape, scale, end))
            } else {
                Ok(Extended::Gamma {
                    p: shape,
                    q: scale,
                    end,
                })
            }
        }
        "TWOGAM" | "TWOGAMpw" => {
            if !(p.len() == 5 || p.len() == 6)
                || p[0] <= 0.0
                || p[1] <= 0.0
                || p[3] <= 0.0
                || p[4] <= 0.0
                || p[2].abs() > 1.0
            {
                return Err(bad());
            }
            let d = p.get(5).copied().unwrap_or(0.0);
            if d < 0.0 {
                return Err(bad());
            }
            let mut a = [p[0], p[1], p[2], p[3], p[4]];
            if name == "TWOGAMpw" {
                (a[0], a[1]) = peak_width(a[0], a[1])?;
                (a[3], a[4]) = peak_width(a[3], a[4])?;
            }
            let end = gamma_end(a[0], a[1], 9.9)?.max(gamma_end(a[3], a[4], 9.9)?);
            if d > 0.0 {
                return table(d, 1, end, |_, t| twogamma(t, a, end));
            }
            // AFNI advances this peak-search coordinate by repeated float
            // addition (not index*DT). Keep that rounding path deliberately.
            let mut t = DT;
            let mut peak = 0.0_f32;
            while t < end {
                peak = peak.max(twogamma(t, a, end).abs());
                t += DT;
            }
            if !peak.is_normal() {
                return Err(bad());
            }
            Ok(Extended::TwoGamma {
                parameters: a,
                end,
                scale: 1.0 / peak,
            })
        }
        "IGFUN" | "IGFUN4" | "IGFUN5" | "UBLOCK" | "UBLOCK4" | "UBLOCK5" | "dmBLOCK"
        | "dmBLOCK4" | "dmBLOCK5" | "dmUBLOCK" | "dmUBLOCK4" | "dmUBLOCK5" => {
            let dynamic = name.starts_with("dm");
            if (dynamic && p.len() > 1) || (!dynamic && !(p.len() == 1 || p.len() == 2)) {
                return Err(bad());
            }
            let duration = if dynamic { 0.0 } else { p[0] };
            let peak = if dynamic {
                p.first().copied().unwrap_or(0.0)
            } else {
                p.get(1).copied().unwrap_or(0.0)
            };
            if !dynamic && (duration <= 0.0 || peak < 0.0) {
                return Err(bad());
            }
            if duration > 999.0 || (dynamic && peak < -999.0) {
                return Err(bad());
            }
            Ok(Extended::Block {
                duration,
                peak,
                power: if name.ends_with('5') { 5 } else { 4 },
                unit: name.contains("UBLOCK"),
                dynamic,
            })
        }
        "CSPLIN" | "CSPLINzero" | "POLY" | "SIN" | "TRIG" => {
            if p.len() != 3 || p[0] >= p[1] {
                return Err(bad());
            }
            // Parse the count from its original text as well: float rounding
            // must not turn a noninteger token into a valid knot count.
            let count = formula
                .rsplit(',')
                .next()
                .and_then(|s| s.strip_suffix(')'))
                .and_then(|s| s.trim().parse::<usize>().ok())
                .ok_or_else(bad)?;
            let minimum = match name {
                "CSPLIN" | "CSPLINzero" => 4,
                "TRIG" => 3,
                _ => 1,
            };
            if count < minimum
                || count > 1 << 24
                || (name == "POLY" && count > 20)
                || !(p[1] - p[0]).is_normal()
            {
                return Err(bad());
            }
            let family = match name {
                "CSPLIN" | "CSPLINzero" => GridFamily::Spline {
                    zero: name.ends_with("zero"),
                },
                "POLY" => GridFamily::Poly,
                "SIN" => GridFamily::Sin,
                _ => GridFamily::Trig,
            };
            Ok(Extended::Grid {
                start: p[0],
                end: p[1],
                count,
                family,
            })
        }
        "SPMG" | "SPMG1" | "SPMG2" | "SPMG3" => {
            let count = match name {
                "SPMG1" => 1,
                "SPMG3" => 3,
                _ => 2,
            };
            if p.is_empty() && !formula.contains('(') {
                return Ok(Extended::Spmg(count));
            }
            if p.len() != 1 || p[0] < 0.0 {
                return Err(bad());
            }
            table(p[0], count, 25.0, spmg)
        }
        "MION" | "MIONN" => {
            let sign = if name == "MIONN" { -1.0 } else { 1.0 };
            if p.is_empty() && !formula.contains('(') {
                return Ok(Extended::Mion(sign));
            }
            if p.len() != 1 || p[0] < 0.0 {
                return Err(bad());
            }
            table(p[0], 1, 60.0, |_, t| sign * mion(t))
        }
        "WAV" => {
            if p.len() > 6 || p.iter().any(|&x| x < 0.0) {
                return Err(bad());
            }
            let mut a = [0.0, 2.0, 4.0, 6.0, 0.2, 2.0];
            for (i, v) in optional.iter().enumerate() {
                if let Some(v) = v {
                    a[i] = *v;
                }
            }
            a[2] = a[2].max(0.1);
            a[3] = a[3].max(0.1);
            a[5] = a[5].max(0.1);
            let end = a[1] + a[2] + a[3] + a[5];
            table(a[0], 1, end, |_, t| wav(t, a))
        }
        _ => Err(bad()),
    }
}

impl Extended {
    pub(crate) fn set_expression_seed(&mut self, seed: u32) -> bool {
        if let Self::Expr { programs, .. } = self {
            for (index, program) in programs.iter().enumerate() {
                // Adjacent components use distinct reproducible streams.
                program.set_seed(seed.wrapping_add(index as u32));
            }
            true
        } else {
            false
        }
    }
    pub(crate) fn components(&self) -> usize {
        match self {
            Self::Grid {
                count,
                family: GridFamily::Spline { zero: true },
                ..
            } => count - 2,
            Self::Grid { count, .. } | Self::Spmg(count) => *count,
            Self::Table { curves, .. } => curves.len(),
            Self::Expr { programs, .. } => programs.len(),
            _ => 1,
        }
    }
    pub(crate) fn duration(&self) -> f32 {
        match self {
            Self::Table { duration, .. } | Self::Block { duration, .. } => *duration,
            _ => 0.0,
        }
    }
    pub(crate) fn requires_duration(&self) -> bool {
        matches!(self, Self::Block { dynamic: true, .. })
    }
    pub(crate) fn support(&self) -> (f32, f32) {
        match self {
            Self::Gamma { end, .. } | Self::TwoGamma { end, .. } => (0.0, *end),
            Self::Block {
                duration, dynamic, ..
            } => (0.0, if *dynamic { 1014.0 } else { duration + 15.0 }),
            Self::Grid { start, end, .. } | Self::Expr { start, end, .. } => (*start, *end),
            Self::Spmg(_) => (0.0, 25.0),
            Self::Mion(_) => (0.0, 60.0),
            // AFNI uses the FIRST component's support for the entire SPMG set.
            Self::Table { curves, .. } => (0.0, DT * curves[0].len() as f32),
        }
    }
    pub(crate) fn sample_duration(&self, t: f32, d: f32) -> Option<f32> {
        if !t.is_finite() || !d.is_finite() || !(0.0..=999.0).contains(&d) {
            return None;
        }
        if let Self::Block {
            peak,
            power,
            unit,
            dynamic: true,
            ..
        } = self
        {
            let v = if d == 0.0 {
                0.0
            } else {
                block(t, d, *peak, *power, *unit)
            };
            v.is_finite().then_some(v)
        } else {
            None
        }
    }
    pub(crate) fn sample(&self, c: usize, t: f32) -> Option<f32> {
        let v = match self {
            Self::Gamma { p, q, end } => gamma(t, *p, *q, *end),
            Self::TwoGamma {
                parameters,
                end,
                scale,
            } => scale * twogamma(t, *parameters, *end),
            Self::Block { dynamic: true, .. } => return None,
            Self::Block {
                duration,
                peak,
                power,
                unit,
                ..
            } => block(t, *duration, *peak, *power, *unit),
            Self::Spmg(_) => spmg(c, t),
            Self::Mion(sign) => sign * mion(t),
            Self::Grid {
                start,
                end,
                count,
                family,
            } => grid(c, t, *start, *end, *count, *family),
            Self::Table { curves, .. } => {
                if t < 0.0 {
                    0.0
                } else {
                    let curve = &curves[c];
                    let x = t / DT;
                    let i = x as usize;
                    let fraction = x - i as f32;
                    if i >= curve.len() {
                        0.0
                    } else {
                        fraction * curve.get(i + 1).copied().unwrap_or(0.0)
                            + (1.0 - fraction) * curve[i]
                    }
                }
            }
            Self::Expr {
                start,
                end,
                programs,
            } => {
                if t < start - 0.0009 || t > end + 0.0009 {
                    0.0
                } else {
                    let x = ((t - start) * (1.0 / (end - start))) as f64;
                    programs[c].evaluate(t as f64, x, 2.0 * x - 1.0)? as f32
                }
            }
        };
        v.is_finite().then_some(v)
    }
}

/// The reference's f32 arithmetic is part of the format's numerical contract.
fn gamma(t: f32, p: f32, q: f32, end: f32) -> f32 {
    if t <= 0.0 || t > end {
        0.0
    } else {
        (((t / (p * q)) as f64).powf(p as f64) * ((p - t / q) as f64).exp()) as f32
    }
}
/// Compute the historical family-specific support cutoff in seconds. Bounds
/// on the grid also prevent a peak scan from running indefinitely when float
/// addition stops advancing at very large coordinates.
fn gamma_end(p: f32, q: f32, multiple: f32) -> Result<f32> {
    let end = p * q + multiple * p.sqrt() * q;
    if !end.is_normal() || end <= 0.0 || end / DT > MAX_GRID as f32 || !(p * q).is_normal() {
        return Err(Error::input(
            "gamma parameters exceed finite supported waveform grid",
        ));
    }
    Ok(end)
}
/// Subtract a delayed gamma undershoot from the primary gamma. The signed
/// ratio is retained; normalization is a separate step in parse_known/table.
fn twogamma(t: f32, p: [f32; 5], end: f32) -> f32 {
    gamma(t, p[0], p[1], end) - p[2] * gamma(t, p[3], p[4], end)
}

/// Both BLOCK conventions use the same integral; only their amplitude rules
/// differ. A negative dmUBLOCK peak means a reference DURATION, not a negative
/// response. A negative dmBLOCK peak leaves its natural amplitude unchanged.
pub(crate) fn block(t: f32, d: f32, peak: f32, power: i32, unit: bool) -> f32 {
    let mut v = block_integral(t, d, power);
    if v > 0.0 && (peak > 0.0 || unit) {
        let reference = if peak > 0.0 {
            d
        } else if peak == 0.0 {
            99.9
        } else {
            -peak
        };
        let tp = (reference as f64 / -(-(reference as f64) / power as f64).exp_m1()) as f32;
        // In AFNI's power-five Maple formula, 1/exp(-tp) overflows
        // for very long blocks. Its normalization guard then skips scaling.
        // Preserve that observable amplitude policy, while our antiderivative
        // still avoids NaN when evaluating the response itself.
        let norm = if power == 5 && !(-(tp as f64)).exp().recip().is_finite() {
            f32::NAN
        } else {
            block_integral(tp, reference, power)
        };
        if norm > 0.0 {
            v *= if peak > 0.0 { peak / norm } else { 1.0 / norm };
        }
    }
    v
}

/// AFNI canonical SPM response, time derivative, and dispersion difference.
/// Component 2 uses the reference 1% perturbation, not a symbolic derivative;
/// the literal coefficients and 25-second endpoint are compatibility choices.
fn spmg(c: usize, x: f32) -> f32 {
    if x <= 0.0 || x >= 25.0 {
        return 0.0;
    }
    let t = x as f64;
    let base = ((-t).exp() * (0.0083333333 * t.powi(5) - 1.274527e-13 * t.powi(15))) as f32;
    match c {
        0 => base,
        1 => {
            ((-t).exp()
                * (0.0083333333 * t.powi(4) * (5.0 - t) - 1.274527e-13 * t.powi(14) * (15.0 - t)))
                as f32
        }
        _ => {
            let perturbed = (0.00869011 * t.powf(4.94057) * (-t / 1.01).exp()
                - 1.274527e-13 * t.powi(15) * (-t).exp()) as f32;
            100.0 * (base - perturbed)
        }
    }
}
/// MION response from the reference exponential mixture, with a 60-second
/// support. MIONN changes only its sign; duration forms integrate this curve.
fn mion(t: f32) -> f32 {
    if t <= 0.0 || t > 60.0 {
        0.0
    } else {
        16.4486
            * (-0.184 / 1.5 * (-t / 1.5).exp()
                + 0.330 / 4.5 * (-t / 4.5).exp()
                + 0.670 / 13.5 * (-t / 13.5).exp())
    }
}
/// Smooth bounded rise used by WAV. The branch endpoints avoid evaluating
/// the tangent expression outside the intended transition interval.
// Retain the source's literal spelling for line-by-line scientific review.
#[allow(clippy::excessive_precision)]
fn transition(x: f32) -> f32 {
    if x <= 0.0 {
        0.0
    } else if x >= 1.0 {
        1.0
    } else {
        let y = std::f64::consts::FRAC_PI_2 * (1.6 * x as f64 - 0.8);
        (0.50212657_f32 as f64 * (y.tan().tanh() + 0.99576486_f32 as f64)) as f32
    }
}
/// Evaluate WAV from onset/rise/fall/undershoot parameters in seconds. The
/// six-element array follows the parser order; zero/default replacement has
/// already happened before this function receives it.
fn wav(t: f32, p: [f32; 6]) -> f32 {
    let rise = p[1];
    let fall = rise + p[2];
    let restore = fall + p[3];
    let end = restore + p[5];
    if t < rise {
        0.0
    } else if t < fall {
        transition((t - rise) / p[2])
    } else if t < restore {
        (1.0 + p[4]) * transition((restore - t) / p[3]) - p[4]
    } else if t < end {
        -p[4] * transition((end - t) / p[5])
    } else {
        0.0
    }
}

/// Reproduce setup_WFUN_function's inclusive, float-accumulated convolution.
/// The inner coordinate is repeatedly SUBTRACTED by DT; replacing it with
/// `(i-j)*DT` or a running sum changes rounding and can change the tail length.
/// Each component is normalized independently to its sampled absolute peak.
fn table<F: Fn(usize, f32) -> f32>(
    duration: f32,
    count: usize,
    end: f32,
    waveform: F,
) -> Result<Extended> {
    if !duration.is_finite() || !(0.0..=999.0).contains(&duration) || !end.is_finite() || end < 0.0
    {
        return Err(Error::input("invalid duration waveform parameters"));
    }
    let lags = 1 + (duration / DT).round_ties_even() as usize;
    let bound = ((duration + end + 2.0) / DT).ceil() as usize;
    if bound > MAX_GRID
        || bound
            .checked_mul(lags)
            .and_then(|n| n.checked_mul(count))
            .is_none_or(|n| n > MAX_CONVOLUTION_TERMS)
    {
        return Err(Error::input(
            "duration waveform exceeds the documented convolution resource limit",
        ));
    }
    let mut curves = Vec::new();
    for c in 0..count {
        let mut curve = Vec::new();
        curve
            .try_reserve(bound)
            .map_err(|_| Error::input("cannot allocate duration waveform"))?;
        let mut peak = 0.001_f32;
        let mut threshold = 0.001_f32;
        let mut nonzero = 0;
        let mut stopped = false;
        for i in 0..bound {
            let mut t = i as f32 * DT;
            let mut value = 0.0_f32;
            for _ in 0..lags {
                if t < 0.0 {
                    break;
                }
                value += waveform(c, t);
                t -= DT;
            }
            if !value.is_finite() {
                return Err(Error::input("nonfinite duration waveform"));
            }
            if value != 0.0 {
                nonzero += 1;
                if value.abs() > peak {
                    peak = value.abs();
                    threshold = 0.0001 * peak;
                }
            } else if i > 99
                && nonzero > 1
                && curve[i - 4..i].iter().all(|v: &f32| v.abs() <= threshold)
            {
                curve.truncate(i - 3);
                stopped = true;
                break;
            }
            curve.push(value);
        }
        if !stopped {
            return Err(Error::input(
                "duration waveform is zero or failed to terminate within its support",
            ));
        }
        let inverse = 1.0 / peak;
        for value in &mut curve {
            *value *= inverse;
        }
        curves.push(curve);
    }
    Ok(Extended::Table { duration, curves })
}

/// Cardinal cubic kernel: one at its own knot, zero at other integer knots,
/// and compact support within two knot spacings. y is distance in knot units.
fn spline(y: f32) -> f32 {
    let y = y.abs();
    if y >= 2.0 {
        0.0
    } else if y >= 1.0 {
        -0.5 * (-4.0 + y * (8.0 + y * (-5.0 + y)))
    } else {
        1.0 + y * y * (1.5 * y - 2.5)
    }
}
/// Map event-relative seconds to a family coordinate. Spline endpoint
/// adjustments preserve AFNI interpolation constraints; zero variants remove
/// endpoint columns rather than modifying the surviving column functions.
fn grid(c: usize, t: f32, start: f32, end: f32, count: usize, family: GridFamily) -> f32 {
    use std::f64::consts::PI;
    match family {
        GridFamily::Spline { zero } => {
            let step = (end - start) / (count - 1) as f32;
            let knot = c + usize::from(zero);
            let y = (t - (start + knot as f32 * step)) / step;
            let lower = if knot == 0 {
                0.0
            } else if knot == 1 {
                -1.0
            } else {
                -2.0
            };
            let upper = if knot == count - 1 {
                0.0
            } else if knot == count - 2 {
                1.0
            } else {
                2.0
            };
            if y < lower - 0.0009 || y > upper + 0.0009 {
                0.0
            } else {
                spline(y)
            }
        }
        GridFamily::Sin => {
            let frequency = (PI * (c + 1) as f64 / (end - start) as f64) as f32;
            if t <= start || t >= end {
                0.0
            } else {
                (frequency * (t - start)).sin()
            }
        }
        GridFamily::Trig => {
            if c == 0 {
                return if t < start - 0.0009 || t > end + 0.0009 {
                    0.0
                } else {
                    1.0
                };
            }
            let frequency = (2.0 * PI * c.div_ceil(2) as f64 / (end - start) as f64) as f32;
            if c.is_multiple_of(2) {
                if t <= start || t >= end {
                    0.0
                } else {
                    (frequency * (t - start)).sin()
                }
            } else if t < start - 0.0009 || t > end + 0.0009 {
                0.0
            } else {
                (frequency * (t - start)).cos()
            }
        }
        GridFamily::Poly => {
            if c == 0 {
                return if t < start - 0.0009 || t > end + 0.0009 {
                    0.0
                } else {
                    1.0
                };
            }
            let x = 2.0 * (t - start) / (end - start) - 1.0;
            if !(-1.000001..=1.000001).contains(&x) {
                0.0
            } else {
                legendre_reference(c, x)
            }
        }
    }
}

/// AFNI's gamma peak/FWHM conversion minimizes the error of two six-step
/// Newton solves. We port that search, rather than substituting a Gaussian
/// width approximation, because the converted f32 parameters enter metadata
/// comparisons and every sample. All state is local (no C static globals).
fn peak_width(peak: f32, width: f32) -> Result<(f32, f32)> {
    if peak <= 0.0 || width <= 0.0 {
        return Err(Error::input("gamma peak and FWHM must be positive"));
    }
    let (p, q) = peak_width_double(peak as f64, width as f64)
        .ok_or_else(|| Error::input("invalid gamma peak/FWHM conversion"))?;
    gamma_end(p as f32, q as f32, 9.9)?;
    Ok((p as f32, q as f32))
}
/// Double-precision peak/FWHM conversion shared with EXPR GAMP/GAMQ. None
/// signals invalid/nonfinite conversion; callers must not substitute defaults.
pub(crate) fn peak_width_double(peak: f64, width: f64) -> Option<(f64, f64)> {
    if !peak.is_finite() || !width.is_finite() || peak <= 0.0 || width <= 0.0 {
        return None;
    }
    let cost = |p: f64| {
        let q = peak / p;
        let guess = 2.3 * p.sqrt() * q;
        let root = |mut x: f64, lower: bool| {
            for _ in 0..6 {
                let dx = (p * (x / (p * q)).ln() + p - x / q + std::f64::consts::LN_2)
                    / (p / x - 1.0 / q);
                x -= dx;
                if lower && x <= 0.0 {
                    x = 0.5 * (x + dx);
                }
                if dx.abs() < 0.000001 * x {
                    break;
                }
            }
            x
        };
        let lo = if peak - 0.5 * guess > 0.0 {
            peak - 0.5 * guess
        } else {
            0.5 * peak
        };
        (root(peak + 0.3 * guess, false) - root(lo, true) - width).abs()
    };
    let initial = (2.3 * peak / width).powi(2);
    let first = minimize(0.5 * initial, 1.5 * initial, &cost);
    let p = minimize(0.8 * first, 1.2 * first, &cost);
    let q = peak / p;
    (p.is_finite() && q.is_finite() && p > 0.0 && q > 0.0).then_some((p, q))
}
/// Bounded one-dimensional search adapted from powell_int.c:minimize_1D.
/// Evaluate 26 equally spaced candidates and refine around the best for five
/// levels. F is a callable cost function; &F borrows it without copying captured
/// parameters. This is a compatibility search, not a general optimizer.
pub(crate) fn minimize<F: Fn(f64) -> f64>(mut lower: f64, mut upper: f64, cost: &F) -> f64 {
    for level in 0..5 {
        let step = (upper - lower) / 25.0;
        let mut best = 0;
        let mut value = f64::INFINITY;
        for i in 0..=25 {
            let v = cost(lower + i as f64 * step);
            if i == 0 || v < value {
                best = i;
                value = v;
            }
        }
        if level == 4 {
            return lower + best as f64 * step;
        }
        if best == 0 {
            upper = lower + 1.5 * step;
        } else if best == 25 {
            lower = upper - 1.5 * step;
        } else {
            lower += (best - 1) as f64 * step * 0.99;
            upper = lower + 1.98 * step;
        }
    }
    unreachable!("the bounded search returns on its fifth level")
}

/// Separate the support interval from expression text, then compile each
/// component once. Programs own their trees; no reference to formula survives.
fn parse_expression(formula: &str) -> Result<Extended> {
    let (interval, text) = formula
        .strip_prefix("EXPR(")
        .and_then(|s| s.split_once(')'))
        .ok_or_else(|| invalid(formula))?;
    let p: Vec<_> = interval.split(',').collect();
    if p.len() != 2 {
        return Err(invalid(formula));
    }
    let start: f32 = p[0].trim().parse().map_err(|_| invalid(formula))?;
    let end: f32 = p[1].trim().parse().map_err(|_| invalid(formula))?;
    if !start.is_finite() || !end.is_finite() || start >= end || !(end - start).is_normal() {
        return Err(invalid(formula));
    }
    let programs = Expression::parse_list(text)?;
    Ok(Extended::Expr {
        start,
        end,
        programs,
    })
}

// The reference uses float explicit polynomials for degrees 0..9, and double
// Horner expressions of float x*x for degrees 10..20. The generated arms below
// retain that evaluation order, including historical coefficient precision.
// Preserve the published coefficient digits, even where a trailing zero
// exceeds what the destination type can distinguish.
#[allow(clippy::excessive_precision)]
fn legendre_reference(n: usize, x: f32) -> f32 {
    let xq = x * x;
    // The pinned ARM64 C build contracts multiply+add in these float Horner
    // expressions. Explicit mul_add gives one rounding on every Rust target;
    // separate operations lose several micro-units near degree-nine roots.
    let coefficients: &[f32] = match n {
        0 => return 1.0,
        1 => return x,
        2 => &[3.0, -1.0],
        3 => &[5.0, -3.0],
        4 => &[35.0, -30.0, 3.0],
        5 => &[63.0, -70.0, 15.0],
        6 => &[231.0, -315.0, 105.0, -5.0],
        7 => &[429.0, -693.0, 315.0, -35.0],
        8 => &[6435.0, -12012.0, 6930.0, -1260.0, 35.0],
        9 => &[12155.0, -25740.0, 18018.0, -4620.0, 315.0],
        _ => &[],
    };
    if !coefficients.is_empty() {
        let mut value = coefficients[0];
        for &coefficient in &coefficients[1..] {
            value = value.mul_add(xq, coefficient);
        }
        let divisor = match n {
            2 | 3 => 2.0,
            4 | 5 => 8.0,
            6 | 7 => 16.0,
            _ => 128.0,
        };
        return value * if n % 2 == 1 { x } else { 1.0 } / divisor;
    }
    match n {
        10 => {
            (-0.24609375e0
                + (0.1353515625e2
                    + (-0.1173046875e3
                        + (0.3519140625e3
                            + (-0.42732421875e3 + 0.18042578125e3 * (xq as f64)) * (xq as f64))
                            * (xq as f64))
                        * (xq as f64))
                    * (xq as f64)) as f32
        }
        11 => {
            ((-0.270703125e1
                + (0.5865234375e2
                    + (-0.3519140625e3
                        + (0.8546484375e3
                            + (-0.90212890625e3 + 0.34444921875e3 * (xq as f64)) * (xq as f64))
                            * (xq as f64))
                        * (xq as f64))
                    * (xq as f64))
                * (x as f64)) as f32
        }
        12 => {
            (0.2255859375e0
                + (-0.17595703125e2
                    + (0.2199462890625e3
                        + (-0.99708984375e3
                            + (0.20297900390625e4
                                + (-0.1894470703125e4 + 0.6601943359375e3 * (xq as f64))
                                    * (xq as f64))
                                * (xq as f64))
                            * (xq as f64))
                        * (xq as f64))
                    * (xq as f64)) as f32
        }
        13 => {
            ((0.29326171875e1
                + (-0.87978515625e2
                    + (0.7478173828125e3
                        + (-0.270638671875e4
                            + (0.47361767578125e4
                                + (-0.3961166015625e4 + 0.12696044921875e4 * (xq as f64))
                                    * (xq as f64))
                                * (xq as f64))
                            * (xq as f64))
                        * (xq as f64))
                    * (xq as f64))
                * (x as f64)) as f32
        }
        14 => {
            (-0.20947265625e0
                + (0.2199462890625e2
                    + (-0.37390869140625e3
                        + (0.2368088378906250e4
                            + (-0.710426513671875e4
                                + (0.1089320654296875e5
                                    + (-0.825242919921875e4
                                        + 0.244852294921875e4 * (xq as f64))
                                        * (xq as f64))
                                    * (xq as f64))
                                * (xq as f64))
                            * (xq as f64))
                        * (xq as f64))
                    * (xq as f64)) as f32
        }
        15 => {
            ((-0.314208984375e1
                + (0.12463623046875e3
                    + (-0.142085302734375e4
                        + (0.7104265136718750e4
                            + (-0.1815534423828125e5
                                + (0.2475728759765625e5
                                    + (-0.1713966064453125e5
                                        + 0.473381103515625e4 * (xq as f64))
                                        * (xq as f64))
                                    * (xq as f64))
                                * (xq as f64))
                            * (xq as f64))
                        * (xq as f64))
                    * (xq as f64))
                * (x as f64)) as f32
        }
        16 => {
            (0.196380615234375e0
                + (-0.26707763671875e2
                    + (0.5920220947265625e3
                        + (-0.4972985595703125e4
                            + (0.2042476226806641e5
                                + (-0.4538836059570312e5
                                    + (0.5570389709472656e5
                                        + (-0.3550358276367188e5
                                            + 0.9171758880615234e4 * (xq as f64))
                                            * (xq as f64))
                                        * (xq as f64))
                                    * (xq as f64))
                                * (xq as f64))
                            * (xq as f64))
                        * (xq as f64))
                    * (xq as f64)) as f32
        }
        17 => {
            ((0.3338470458984375e1
                + (-0.1691491699218750e3
                    + (0.2486492797851562e4
                        + (-0.1633980981445312e5
                            + (0.5673545074462891e5
                                + (-0.1114077941894531e6
                                    + (0.1242625396728516e6
                                        + (-0.7337407104492188e5
                                            + 0.1780400253295898e5 * (xq as f64))
                                            * (xq as f64))
                                        * (xq as f64))
                                    * (xq as f64))
                                * (xq as f64))
                            * (xq as f64))
                        * (xq as f64))
                    * (xq as f64))
                * (x as f64)) as f32
        }
        18 => {
            (-0.1854705810546875e0
                + (0.3171546936035156e2
                    + (-0.8880331420898438e3
                        + (0.9531555725097656e4
                            + (-0.5106190567016602e5
                                + (0.1531857170104980e6
                                    + (-0.2692355026245117e6
                                        + (0.2751527664184570e6
                                            + (-0.1513340215301514e6
                                                + 0.3461889381408691e5 * (xq as f64))
                                                * (xq as f64))
                                            * (xq as f64))
                                        * (xq as f64))
                                    * (xq as f64))
                                * (xq as f64))
                            * (xq as f64))
                        * (xq as f64))
                    * (xq as f64)) as f32
        }
        19 => {
            ((-0.3523941040039062e1
                + (0.2220082855224609e3
                    + (-0.4084952453613281e4
                        + (0.3404127044677734e5
                            + (-0.1531857170104980e6
                                + (0.4038532539367676e6
                                    + (-0.6420231216430664e6
                                        + (0.6053360861206055e6
                                            + (-0.3115700443267822e6
                                                + 0.6741574058532715e5 * (xq as f64))
                                                * (xq as f64))
                                            * (xq as f64))
                                        * (xq as f64))
                                    * (xq as f64))
                                * (xq as f64))
                            * (xq as f64))
                        * (xq as f64))
                    * (xq as f64))
                * (x as f64)) as f32
        }
        20 => {
            (0.1761970520019531e0
                + (-0.3700138092041016e2
                    + (0.1276547641754150e4
                        + (-0.1702063522338867e5
                            + (0.1148892877578735e6
                                + (-0.4442385793304443e6
                                    + (0.1043287572669983e7
                                        + (-0.1513340215301514e7
                                            + (0.1324172688388824e7
                                                + (-0.6404495355606079e6
                                                    + 0.1314606941413879e6 * (xq as f64))
                                                    * (xq as f64))
                                                * (xq as f64))
                                            * (xq as f64))
                                        * (xq as f64))
                                    * (xq as f64))
                                * (xq as f64))
                            * (xq as f64))
                        * (xq as f64))
                    * (xq as f64)) as f32
        }
        _ => unreachable!("POLY count is validated before sampling"),
    }
}
