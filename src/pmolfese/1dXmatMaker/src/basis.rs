// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! AFNI response-function catalog, in seconds.
//!
//! These formulas and compatibility conventions follow 3dDeconvolve.c
//! (basis_gam, basis_tent, basis_spmg1/2, basis_block*_OLD, basis_parser),
//! revision 2664fdd5de7425f658c6d901765fbb4c291693fa. AFNI's source credits
//! B. Douglas Ward and the Medical College of Wisconsin, copyright 1998-2003,
//! released under CC BY 4.0 (<https://creativecommons.org/licenses/by/4.0/>).
//! This Rust adaptation changes storage/validation and expresses the BLOCK
//! integral as an exponential polynomial rather than Maple-generated terms.
//!
//! Float intermediates here are intentional: AFNI samples basis functions in
//! f32. We widen completed samples into the f64 design, not the reverse. Read
//! parse(), then sample(); the additional families and duration tables live
//! in basis_extended.rs so the elementary formulas below stay readable.
//!
//! `Kind` holds the family-specific parameters; Basis also owns the exact
//! formula string for provenance. The boxed Extended variant stores a pointer
//! in this enum instead of making every simple GAM carry its largest payload.
//! Box is an owning allocation, not a borrowed C pointer: Rust frees it with
//! Basis. The numerical routines borrow it during sampling without copying it.

use crate::error::{Error, Result};

/// A parsed response family; private fields protect parameter bounds.
/// EXPR random calls advance a private stream; other sampling is immutable.
#[derive(Debug)]
pub struct Basis {
    formula: String,
    kind: Kind,
}

#[derive(Debug)]
enum Kind {
    Sampled(crate::sampled::Curve),
    Extended(Box<crate::basis_extended::Extended>),
    Gamma {
        shape: f32,
        scale: f32,
        end: f32,
    },
    Block {
        duration: f32,
        peak: f32,
        power: i32,
    },
    Tent {
        start: f32,
        end: f32,
        knots: usize,
        zero: bool,
    },
    Spmg {
        components: usize,
    },
}

impl Basis {
    /// Own an already prepared sampled HRF; its fixed duration must not be
    /// reapplied by timing assembly. The descriptive formula is provenance.
    pub fn sampled(curve: crate::sampled::Curve, formula: String) -> Result<Self> {
        if curve.end() > f32::MAX as f64
            || formula.is_empty()
            || formula.chars().any(char::is_control)
        {
            return Err(Error::input("invalid sampled basis support/provenance"));
        }
        Ok(Self {
            formula,
            kind: Kind::Sampled(curve),
        })
    }

    /// Parse the documented AFNI basis catalog, including fixed duration and
    /// EXPR forms. Unknown names and malformed parameters never fall back.
    /// Parameters must remain finite in AFNI's float representation. Negative
    /// BLOCK peaks are rejected rather than silently treated as zero.
    pub fn parse(formula: &str) -> Result<Self> {
        if let Some(result) = crate::basis_extended::parse(formula) {
            let mut extended = result?;
            extended.set_expression_seed(0x1234_abcd);
            return Ok(Self {
                formula: formula.to_owned(),
                kind: Kind::Extended(Box::new(extended)),
            });
        }
        let bad = || {
            Error::input(format!(
                "unsupported or invalid basis '{formula}'; see -help"
            ))
        };
        let (name, parameters) = if let Some((name, rest)) = formula.split_once('(') {
            let body = rest.strip_suffix(')').ok_or_else(bad)?;
            let mut numbers = Vec::new();
            for part in body.split(',') {
                let number: f32 = part.trim().parse().map_err(|_| bad())?;
                if !number.is_finite() {
                    return Err(bad());
                }
                numbers.push(number);
            }
            (name, numbers)
        } else {
            (formula, Vec::new())
        };
        let kind = match (name, parameters.as_slice()) {
            ("GAM", []) if formula == "GAM" => Kind::Gamma {
                shape: 8.6,
                scale: 0.547,
                end: 11.1,
            },
            ("GAM", &[shape, scale]) if shape > 0.0 && scale > 0.0 => {
                let end = shape * scale + 9.9 * shape.sqrt() * scale;
                if !end.is_finite() || end <= 0.0 || !(shape * scale).is_normal() {
                    return Err(bad());
                }
                Kind::Gamma { shape, scale, end }
            }
            ("BLOCK" | "BLOCK4" | "BLOCK5", [duration, ..])
                if (parameters.len() == 1 || parameters.len() == 2) && *duration > 0.0 =>
            {
                let peak = parameters.get(1).copied().unwrap_or(0.0);
                if peak < 0.0 || !(*duration + 15.0).is_finite() {
                    return Err(bad());
                }
                let power = if name == "BLOCK5" { 5 } else { 4 };
                let tp = (*duration as f64 / -(-(*duration as f64) / power as f64).exp_m1()) as f32;
                if !tp.is_finite() || block_integral(tp, *duration, power) <= 0.0 {
                    return Err(bad());
                }
                Kind::Block {
                    duration: *duration,
                    peak,
                    power,
                }
            }
            ("TENT" | "TENTzero", &[start, end, n])
                if start < end && n.fract() == 0.0 && n >= 2.0 && n < i32::MAX as f32 =>
            {
                // The knot count is an integer token, not a float that may
                // round 2.00000001 down to 2. Require exact f32 integer range
                // too, since AFNI's knot-coordinate arithmetic is float.
                let knots = formula
                    .rsplit(',')
                    .next()
                    .and_then(|s| s.strip_suffix(')'))
                    .and_then(|s| s.trim().parse::<usize>().ok())
                    .ok_or_else(bad)?;
                if knots > (1 << 24) || knots != n as usize {
                    return Err(bad());
                }
                let zero = name == "TENTzero";
                let step = (end - start) / (n - 1.0);
                if (zero && knots < 3)
                    || !step.is_normal()
                    || !(start - 0.00111 * step).is_finite()
                    || !(end + 0.00111 * step).is_finite()
                {
                    return Err(bad());
                }
                Kind::Tent {
                    start,
                    end,
                    knots,
                    zero,
                }
            }
            ("SPMG1" | "SPMG2", []) if !formula.contains('(') => Kind::Spmg {
                components: if name == "SPMG1" { 1 } else { 2 },
            },
            _ => return Err(bad()),
        };
        Ok(Self {
            formula: formula.to_owned(),
            kind,
        })
    }

    /// Reset EXPR random streams; return false for a non-EXPR family.
    /// The default base seed is 0x1234abcd. Component c uses seed+c modulo
    /// 2^32. Subsequent sampling consumes values in call order. This is a
    /// reproducible local policy, not AFNI's clock/PID-dependent sequence.
    /// `&mut self` requires exclusive access when changing the seed.
    pub fn set_expression_seed(&mut self, seed: u32) -> bool {
        match &mut self.kind {
            Kind::Extended(basis) => basis.set_expression_seed(seed),
            _ => false,
        }
    }

    /// Duration-modulated families need a duration from each event. Their
    /// waveform API is available in M3; married timing/AM assembly is M4.
    pub fn requires_duration(&self) -> bool {
        matches!(&self.kind, Kind::Extended(b) if b.requires_duration())
    }
    /// Evaluate a dmBLOCK/dmUBLOCK waveform with an explicit event duration.
    /// Durations are seconds in `[0,999]`; zero produces an all-zero response.
    /// Other families reject this call instead of ignoring the supplied value.
    pub fn sample_with_duration(&self, seconds: f64, duration: f64) -> Option<f64> {
        match &self.kind {
            Kind::Extended(b) if seconds.is_finite() && duration.is_finite() => b
                .sample_duration(seconds as f32, duration as f32)
                .map(f64::from),
            _ => None,
        }
    }
    /// Preserve the user's exact supported formula for BasisFormula metadata.
    pub fn formula(&self) -> &str {
        &self.formula
    }
    /// Number of coefficients per event (or per condition for ordinary timing).
    pub fn components(&self) -> usize {
        match self.kind {
            Kind::Extended(ref basis) => basis.components(),
            Kind::Tent { knots, zero, .. } => knots - if zero { 2 } else { 0 },
            Kind::Spmg { components } => components,
            _ => 1,
        }
    }
    /// Fixed event duration; zero denotes an impulse, not an HRF support length.
    pub fn duration(&self) -> f64 {
        match self.kind {
            Kind::Extended(ref basis) => basis.duration() as f64,
            Kind::Block { duration, .. } => duration as f64,
            Kind::Sampled(ref curve) => curve.duration(),
            _ => 0.0,
        }
    }
    /// Closed model support used by the sampling loop; basis endpoint rules
    /// may additionally return zero at either endpoint (e.g. SPMG at 25 s).
    pub(crate) fn support(&self) -> (f32, f32) {
        match self.kind {
            Kind::Extended(ref basis) => basis.support(),
            Kind::Gamma { end, .. } => (0.0, end),
            Kind::Block { duration, .. } => (0.0, duration + 15.0),
            Kind::Tent { start, end, .. } => (start, end),
            Kind::Spmg { .. } => (0.0, 25.0),
            Kind::Sampled(ref curve) => (0.0, curve.end() as f32),
        }
    }
    /// Evaluate one component at seconds relative to onset. None means invalid
    /// component/time or a nonfinite numerical result, never a substituted model.
    pub fn sample(&self, component: usize, seconds: f64) -> Option<f64> {
        if component >= self.components() || !seconds.is_finite() {
            return None;
        }
        let x = seconds as f32;
        if !x.is_finite() {
            return None;
        }
        let value = match self.kind {
            Kind::Extended(ref basis) => basis.sample(component, x)?,
            Kind::Sampled(ref curve) => curve.sample(x as f64)? as f32,
            Kind::Gamma { shape, scale, end } => {
                if x <= 0.0 || x > end {
                    0.0
                } else {
                    // The ratio/exponent subtraction are float in the C source;
                    // pow and exp then use double before returning float.
                    let power = ((x / (shape * scale)) as f64).powf(shape as f64);
                    let decay = ((shape - x / scale) as f64).exp();
                    (power * decay) as f32
                }
            }
            Kind::Block {
                duration,
                peak,
                power,
            } => {
                // The shared Rust helper keeps aliases and dm/fixed forms on
                // one amplitude policy. This is not an afni-core/FFI call.
                crate::basis_extended::block(x, duration, peak, power, false)
            }
            Kind::Tent {
                start,
                end,
                knots,
                zero,
            } => {
                let step = (end - start) / (knots - 1) as f32;
                let knot = component + usize::from(zero);
                let mid = if knot == knots - 1 {
                    end
                } else {
                    start + knot as f32 * step
                };
                let left = if knot == 0 {
                    start - 0.00111 * step
                } else {
                    start + (knot - 1) as f32 * step
                };
                let right = if knot == knots - 1 {
                    end + 0.00111 * step
                } else {
                    start + (knot + 1) as f32 * step
                };
                if x <= left || x >= right {
                    0.0
                } else if x <= mid {
                    (x - left) / (mid - left)
                } else {
                    (right - x) / (right - mid)
                }
            }
            Kind::Spmg { .. } => {
                if x <= 0.0 || x >= 25.0 {
                    0.0
                } else {
                    let t = x as f64;
                    // The derivative is a distinct, unnormalized coefficient;
                    // do not normalize either curve or orthogonalize the pair.
                    let y = if component == 0 {
                        0.0083333333 * t.powi(5) - 1.274527e-13 * t.powi(15)
                    } else {
                        0.0083333333 * t.powi(4) * (5.0 - t)
                            - 1.274527e-13 * t.powi(14) * (15.0 - t)
                    };
                    ((-t).exp() * y) as f32
                }
            }
        };
        value.is_finite().then_some(value as f64)
    }
}

/// Integrate the peak-one gamma kernel over a box of duration d.
/// For integer n, an antiderivative of t^n exp(-t) is
/// -n! exp(-t) sum(k=0..n, t^k/k!). Subtract its values at max(0,t-d), t.
/// Keep AFNI's rounded exp(4)/exp(5) constants and t>=d+15 cutoff: replacing
/// these with exact constants changes the historical natural-amplitude model.
pub(crate) fn block_integral(seconds: f32, duration: f32, power: i32) -> f32 {
    if seconds <= 0.0 || seconds >= duration + 15.0 {
        return 0.0;
    }
    let t = seconds as f64;
    let low = (t - duration as f64).max(0.0);
    let polynomial = |x: f64| {
        let mut term = 1.0;
        let mut sum = 1.0;
        for k in 1..=power {
            term *= x / k as f64;
            sum += term;
        }
        sum
    };
    let (factorial, scale) = if power == 4 {
        (24.0, 256.0)
    } else {
        (120.0, 3125.0)
    };
    let (upper_exp, lower_exp) = if power == 4 {
        (
            (4.0 - t).exp(),
            if t < duration as f64 {
                54.5982
            } else {
                (4.0 - low).exp()
            },
        )
    } else {
        (148.413 * (-t).exp(), 148.413 * (-low).exp())
    };
    (factorial / scale * (lower_exp * polynomial(low) - upper_exp * polynomial(t))) as f32
}
