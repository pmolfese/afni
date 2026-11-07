// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Scalar numerical support for EXPR, with no native-library dependency.
//!
//! AFNI parser_int.c/mri_stats.c supply the public domain/error conventions
//! (MCW 1994–2000, CC BY 4.0). Probability formulas below use defining series,
//! continued fractions, and bracketed inversion, not AFNI's Fortran/C runtime.
//! A finite iteration budget reports nonconvergence as NaN; the expression
//! evaluator turns that into a contextual error before publishing a matrix.
//!
//! Why separate this module? Parsing decides *what* to evaluate. These helpers
//! decide *how* to evaluate the mathematical functions accurately. Keeping the
//! two apart makes reference comparisons and independent identities useful.

use std::f64::consts::{PI, SQRT_2};
const EPS: f64 = 2.0e-14;
const TINY: f64 = 1.0e-300;
const BIG: f64 = 9.99e37;

/// Log Gamma for positive arguments. Shift to the Stirling region using
/// Gamma(x+1)=x Gamma(x), then retain Bernoulli corrections through x^-15.
fn log_gamma(mut x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    let mut correction = 0.0;
    while x < 12.0 {
        correction -= x.ln();
        x += 1.0;
    }
    let z = 1.0 / (x * x);
    correction + (x - 0.5) * x.ln() - x
        + 0.5 * (2.0 * PI).ln()
        + (1.0 / 12.0
            + z * (-1.0 / 360.0
                + z * (1.0 / 1260.0
                    + z * (-1.0 / 1680.0
                        + z * (1.0 / 1188.0
                            + z * (-691.0 / 360360.0
                                + z * (1.0 / 156.0 - z * 3617.0 / 122400.0)))))))
            / x
}
/// Extend positive Gamma to negative nonintegers with Euler reflection.
/// Poles return NaN so the expression evaluator rejects an undefined sample.
pub(crate) fn gamma(x: f64) -> f64 {
    if x > 0.0 {
        log_gamma(x).exp()
    } else if x == x.trunc() {
        f64::NAN
    } else {
        PI / ((PI * x).sin() * log_gamma(1.0 - x).exp())
    }
}

/// Return BOTH regularized gamma tails. Computing the small tail directly
/// avoids losing it to subtraction from 1 in significance calculations.
fn gamma_pq(a: f64, x: f64) -> (f64, f64) {
    if a <= 0.0 || x < 0.0 {
        return (f64::NAN, f64::NAN);
    }
    if x == 0.0 {
        return (0.0, 1.0);
    }
    let factor = (a * x.ln() - x - log_gamma(a)).exp();
    if x < a + 1.0 {
        let mut term = 1.0 / a;
        let mut sum = term;
        for n in 1..=100_000 {
            term *= x / (a + n as f64);
            sum += term;
            if term.abs() < sum.abs() * EPS {
                let p = (sum * factor).clamp(0.0, 1.0);
                return (p, 1.0 - p);
            }
        }
    } else {
        // Modified Lentz evaluation rescales numerator and denominator rather
        // than forming enormous continued-fraction polynomials explicitly.
        let mut b = x + 1.0 - a;
        let mut c = 1.0 / TINY;
        let mut d = 1.0 / b;
        let mut h = d;
        for n in 1..=100_000 {
            let n = n as f64;
            let an = -n * (n - a);
            b += 2.0;
            d = an * d + b;
            if d.abs() < TINY {
                d = TINY;
            }
            c = b + an / c;
            if c.abs() < TINY {
                c = TINY;
            }
            d = 1.0 / d;
            let step = d * c;
            h *= step;
            if (step - 1.0).abs() < EPS {
                let q = (factor * h).clamp(0.0, 1.0);
                return (1.0 - q, q);
            }
        }
    }
    (f64::NAN, f64::NAN)
}
/// Evaluate the incomplete-beta continued fraction using modified Lentz.
/// Each iteration has TWO substeps; convergence is tested only after both.
/// TINY prevents a zero denominator, while EPS is a relative stopping bound.
fn beta_fraction(a: f64, b: f64, x: f64) -> f64 {
    let mut c = 1.0;
    let mut d = 1.0 - (a + b) * x / (a + 1.0);
    if d.abs() < TINY {
        d = TINY;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..=100_000 {
        let m = m as f64;
        for (substep, aa) in [
            m * (b - m) * x / ((a + 2.0 * m - 1.0) * (a + 2.0 * m)),
            -(a + m) * (a + b + m) * x / ((a + 2.0 * m) * (a + 2.0 * m + 1.0)),
        ]
        .into_iter()
        .enumerate()
        {
            d = 1.0 + aa * d;
            if d.abs() < TINY {
                d = TINY;
            }
            c = 1.0 + aa / c;
            if c.abs() < TINY {
                c = TINY;
            }
            d = 1.0 / d;
            let step = d * c;
            h *= step;
            if substep == 1 && (step - 1.0).abs() < EPS {
                return h;
            }
        }
    }
    f64::NAN
}
/// Regularized incomplete beta and its complementary tail. Reflect about
/// x=1 when needed so the smaller tail is computed directly, not subtracted.
fn beta_pq(x: f64, a: f64, b: f64) -> (f64, f64) {
    if a <= 0.0 || b <= 0.0 {
        return (f64::NAN, f64::NAN);
    }
    if x <= 0.0 {
        return (0.0, 1.0);
    }
    if x >= 1.0 {
        return (1.0, 0.0);
    }
    let factor =
        (log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * x.ln() + b * (-x).ln_1p()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        let p = (factor * beta_fraction(a, b, x) / a).clamp(0.0, 1.0);
        (p, 1.0 - p)
    } else {
        let q = (factor * beta_fraction(b, a, 1.0 - x) / b).clamp(0.0, 1.0);
        (1.0 - q, q)
    }
}
/// erf(x)=sign(x)*P(1/2,x^2), tying error functions to the gamma integral.
pub(crate) fn erf(x: f64) -> f64 {
    x.signum() * gamma_pq(0.5, x * x).0
}
/// Compute the complementary error function directly in the small tail.
pub(crate) fn erfc(x: f64) -> f64 {
    let q = gamma_pq(0.5, x * x).1;
    if x < 0.0 { 2.0 - q } else { q }
}
/// Upper-tail probability `Pr[N(0,1)>x]`, not the ordinary lower-tail CDF.
pub(crate) fn normal_q(x: f64) -> f64 {
    0.5 * erfc(x / SQRT_2)
}
/// Invert the upper normal tail with the reference initial approximation
/// and three Newton refinements. The historical endpoint cap is 13 sigma.
pub(crate) fn normal_q_inverse(p: f64) -> f64 {
    // Preserve parser.c's endpoint convention and its three Newton updates.
    let tail = if p > 0.5 { 1.0 - p } else { p };
    let mut x = 13.0;
    if tail > 0.0 {
        let t = (-2.0 * tail.ln()).sqrt();
        x = t
            - ((t * 0.010328 + 0.802853) * t + 2.525517)
                / (((t * 0.001308 + 0.189269) * t + 1.432788) * t + 1.0);
        for _ in 0..3 {
            x += (normal_q(x) - tail) / ((-0.5 * x * x).exp() / 2.506628274631);
        }
    }
    if p > 0.5 { -x } else { x }
}

/// J0/J1 and Y0/Y1 from their convergent small-argument series and large-
/// argument asymptotic expansions. Stop the asymptotic series before terms
/// grow: it is an asymptotic expansion, not a convergent infinite sum.
pub(crate) fn bessel(input: f64, n: usize, second: bool) -> f64 {
    if second && input <= 0.0 {
        return 0.0;
    }
    let x = input.abs();
    if x == 0.0 {
        return if n == 0 { 1.0 } else { 0.0 };
    }
    let value = if x < 12.0 {
        let z = -0.25 * x * x;
        let mut term = 1.0;
        let mut j0 = 1.0;
        let mut sum_h = 0.0;
        let mut derivative_h = 0.0;
        let mut harmonic = 0.0;
        for k in 1..=200 {
            let k = k as f64;
            term *= z / (k * k);
            harmonic += 1.0 / k;
            j0 += term;
            sum_h -= harmonic * term;
            derivative_h -= harmonic * term * 2.0 * k / x;
            if term.abs() < EPS * j0.abs() {
                break;
            }
        }
        let mut term = x / 2.0;
        let mut j1 = term;
        for k in 1..=200 {
            let k = k as f64;
            term *= z / (k * (k + 1.0));
            j1 += term;
            if term.abs() < EPS * j1.abs() {
                break;
            }
        }
        if !second {
            if n == 0 { j0 } else { j1 }
        } else {
            let log = (x / 2.0).ln() + std::f64::consts::EULER_GAMMA;
            if n == 0 {
                2.0 / PI * (log * j0 + sum_h)
            } else {
                2.0 / PI * (log * j1 - j0 / x - derivative_h)
            }
        }
    } else {
        let mu = 4.0 * (n * n) as f64;
        let mut term = 1.0;
        let mut p = 1.0;
        let mut q = 0.0;
        for k in 1..=100 {
            let next = term * (mu - ((2 * k - 1) * (2 * k - 1)) as f64) / (8.0 * x * k as f64);
            if next.abs() > term.abs() {
                break;
            }
            term = next;
            match k % 4 {
                0 => p += term,
                1 => q += term,
                2 => p -= term,
                _ => q -= term,
            }
            if term.abs() < EPS {
                break;
            }
        }
        let phase = x - (2 * n + 1) as f64 * PI / 4.0;
        (2.0 / (PI * x)).sqrt()
            * if second {
                p * phase.sin() + q * phase.cos()
            } else {
                p * phase.cos() - q * phase.sin()
            }
    };
    if !second && n == 1 && input < 0.0 {
        -value
    } else {
        value
    }
}

/// EXPR PLEG uses a double three-term recurrence. This differs from the
/// POLY basis compatibility evaluator, whose input and low orders are float.
pub(crate) fn legendre(n: i32, x: f64) -> f64 {
    if n < 1 {
        return 1.0;
    }
    if n > 100_000 {
        return f64::NAN;
    }
    let mut a = 1.0;
    let mut b = x;
    for k in 2..=n {
        let next = ((2 * k - 1) as f64 * x * b - (k - 1) as f64 * a) / k as f64;
        a = b;
        b = next;
    }
    b
}
/// EXPR HRFBK uses double inputs and normalizes to unit peak, unlike a bare
/// natural-amplitude BLOCK. Preserve the historical rounded exp constants.
/// EXPR HRFBK4/5: normalized integral of the power-four/five gamma kernel.
/// This is the double expression function, not the float BLOCK basis wrapper.
pub(crate) fn hrf(t: f64, d: f64, power: i32) -> f64 {
    let integral = |t: f64| {
        if t <= 0.0 || t >= d + 15.0 {
            return 0.0;
        }
        let low = (t - d).max(0.0);
        let poly = |x: f64| {
            let mut term = 1.0;
            let mut sum = 1.0;
            for k in 1..=power {
                term *= x / k as f64;
                sum += term;
            }
            sum
        };
        let (fac, div, e) = if power == 4 {
            (24.0, 256.0, 54.5982)
        } else {
            (120.0, 3125.0, 148.413)
        };
        let upper = if power == 4 {
            (4.0 - t).exp()
        } else {
            e * (-t).exp()
        };
        let lower = if power == 4 && t >= d {
            (4.0 - low).exp()
        } else {
            e * (-low).exp()
        };
        fac / div * (lower * poly(low) - upper * poly(t))
    };
    let v = integral(t);
    if v > 0.0 && d > 0.0 {
        v / integral(d / -(-d / power as f64).exp_m1())
    } else {
        v
    }
}
/// Compact rhombic-dodecahedron autocorrelation polynomial from parser_int.c.
/// Sorting absolute coordinates exploits symmetry before selecting each piece.
pub(crate) fn rhdd(a: &[f64]) -> f64 {
    let mut v = [a[0].abs(), a[1].abs(), a[2].abs()];
    v.sort_by(|a, b| b.total_cmp(a));
    let (x, y, z) = (2.0 * v[0], 2.0 * v[1], 2.0 * v[2]);
    let t = x + y - 4.0;
    if v[0] >= 2.0 || t >= 0.0 {
        return 0.0;
    }
    let (xz, yz, xy) = (x + z - 2.0, y + z - 2.0, t + 2.0);
    let pa = 0.00260416666667
        * t.powi(3)
        * (-3.0 * x * y - 5.0 * z * z + 2.0 * (x + y) + 20.0 * z + x * x + y * y - 24.0);
    let pb1 = 0.00520833333333
        * xz.powi(3)
        * (x * x - 9.0 * x - 3.0 * x * z + 10.0 * y - 5.0 * y * y + 14.0 + 11.0 * z + z * z);
    let pb2 = 0.00520833333333
        * yz.powi(3)
        * (46.0 - 30.0 * x - z - y + 3.0 * z * y + 5.0 * x * x - y * y - z * z);
    if xy <= 0.0 {
        pa + pb1 + pb2
            - 0.01041666666667
                * xy.powi(3)
                * (x * x + x - 3.0 * x * y - 5.0 * z * z + y * y + y - 6.0)
    } else if xz <= 0.0 {
        pa + pb1 + pb2
    } else if yz <= 0.0 {
        if x - z >= 2.0 {
            0.00260416666667
                * t.powi(3)
                * (-x * x + 8.0 * x + 3.0 * x * y - y * y + 5.0 * z * z - 16.0 - 12.0 * y)
        } else {
            pa + pb2
        }
    } else {
        pa
    }
}
/// Full width at a requested height of the Gaussian/exponential ACF mixture.
/// The bounded minimizer matches the three reference search refinements.
pub(crate) fn acf_width(a: &[f64]) -> f64 {
    let (mix, b, c, level) = (a[0], a[1], a[2], a[3]);
    if !(0.0..=1.0).contains(&mix) || b <= 0.0 || c <= 0.0 || level <= 0.0 || level >= 1.0 {
        return 0.0;
    }
    let cost = |r: f64| {
        (mix * (-0.5 * r * r / (b * b)).exp() + (1.0 - mix) * (-r / c).exp() - level).abs()
    };
    let high = if level > 0.1 {
        2.5 * b + 2.0 * c
    } else {
        (-5.0 * level.ln()).sqrt() * b + 2.0 * c
    };
    let low = if level > 0.1 {
        0.0666 * high
    } else {
        0.0333 * high
    };
    let r = crate::basis_extended::minimize(low, high, &cost);
    let r = crate::basis_extended::minimize(0.5 * r, 2.0 * r, &cost);
    2.0 * crate::basis_extended::minimize(0.8 * r, 1.25 * r, &cost)
}

/// Invert a monotone CDF using a bracket, not a derivative. Bracketing remains
/// reliable at flat tails and for fractional binomial/Poisson arguments used
/// by AFNI's continuous inverse conventions.
fn inverse<F: Fn(f64) -> f64>(p: f64, mut lo: f64, mut hi: f64, cdf: F) -> f64 {
    if p <= 0.0 {
        return lo;
    }
    if p >= 1.0 {
        return hi;
    }
    for _ in 0..128 {
        if cdf(hi) >= p {
            break;
        }
        hi = 2.0 * hi.abs() + 1.0;
    }
    for _ in 0..128 {
        if cdf(lo) <= p {
            break;
        }
        lo = -2.0 * lo.abs() - 1.0;
    }
    // Expanding a bracket is bounded work, not permission to use an
    // unbracketed answer. Fail explicitly if either endpoint is still wrong.
    if !cdf(lo).is_finite() || !cdf(hi).is_finite() || cdf(lo) > p || cdf(hi) < p {
        return f64::NAN;
    }
    for _ in 0..160 {
        let mid = lo + (hi - lo) / 2.0;
        let value = cdf(mid);
        if !value.is_finite() {
            return f64::NAN;
        }
        if value < p {
            lo = mid;
        } else {
            hi = mid;
        }
        if (hi - lo).abs() <= EPS * (1.0 + mid.abs()) {
            return (lo + hi) / 2.0;
        }
    }
    f64::NAN
}
/// Student-t lower/upper tails through the incomplete beta identity.
/// Reflection around zero preserves the small tail on either side.
fn student_pq(t: f64, df: f64) -> (f64, f64) {
    let tail = 0.5 * beta_pq(df / (df + t * t), df / 2.0, 0.5).0;
    if t < 0.0 {
        (tail, 1.0 - tail)
    } else {
        (1.0 - tail, tail)
    }
}

/// Poisson mixtures define noncentral chi-square/F CDFs. Starting near the
/// Poisson mode avoids exp(-lambda) underflow at large noncentrality.
fn mixture<F: Fn(usize) -> f64>(lambda: f64, value: F) -> f64 {
    if lambda < 0.0 || lambda > 100_000.0 {
        return f64::NAN;
    }
    if lambda == 0.0 {
        return value(0);
    }
    let center = lambda.floor() as usize;
    let weight = (-lambda + center as f64 * lambda.ln() - log_gamma(center as f64 + 1.0)).exp();
    let mut sum = weight * value(center);
    let mut w = weight;
    for j in (0..center).rev() {
        w *= (j + 1) as f64 / lambda;
        sum += w * value(j);
        if w < EPS && (j as f64) < lambda / 2.0 {
            break;
        }
    }
    w = weight;
    for j in center + 1..center + 100_000 {
        w *= lambda / j as f64;
        sum += w * value(j);
        if w < EPS {
            return sum;
        }
    }
    f64::NAN
}
/// Noncentral t uses its beta-series definition. Two Poisson-like weights
/// represent the symmetric and signed parts; negative t follows symmetry.
fn noncentral_t(t: f64, df: f64, delta: f64) -> f64 {
    if df <= 0.0 {
        return 0.0;
    }
    if delta.abs() < 1e-8 {
        return student_pq(t, df).0;
    }
    if t < 0.0 {
        return 1.0 - noncentral_t(-t, df, -delta);
    }
    if t == 0.0 {
        return normal_q(delta);
    }
    let y = t * t / (df + t * t);
    let lambda = delta * delta / 2.0;
    let even = mixture(lambda, |j| beta_pq(y, j as f64 + 0.5, df / 2.0).0);
    let odd = mixture(lambda, |j| {
        (log_gamma(j as f64 + 1.0) - log_gamma(j as f64 + 1.5)).exp()
            * beta_pq(y, j as f64 + 1.0, df / 2.0).0
    });
    (normal_q(delta) + 0.5 * (even + delta / SQRT_2 * odd)).clamp(0.0, 1.0)
}

/// Distinguish invalid distribution parameters (AFNI's documented sentinel)
/// from failure of a numerical algorithm on valid parameters (an error).
/// Without this check, NaN from a failed series could silently become CDF=0.
fn valid_distribution(code: i32, a: f64, b: f64, c: f64) -> bool {
    match code {
        2 | 3 | 6 | 10 | 19 => a > 0.0,
        4 | 7 | 9 | 20 => a > 0.0 && b > 0.0,
        8 => a > 0.0 && (0.0..=1.0).contains(&b),
        11 | 14 | 15 | 21 => b > 0.0,
        12 => a > 0.0 && b > 0.0 && c >= 0.0,
        13 => a > 0.0 && b >= 0.0,
        16 => b > a,
        17 => a > 0.0,
        18 => b > 0.0 && c > 0.0,
        5 | 22..=24 => true,
        _ => false,
    }
}

/// NIfTI statistical intent codes are a numeric API used by STAT2CDF/CDF2STAT.
/// Parameters follow nifticdf.c: notably gamma's second parameter is a RATE.
fn cdf(x: f64, code: i32, a: f64, b: f64, c: f64) -> f64 {
    if !valid_distribution(code, a, b, c) {
        return 0.0;
    }
    let invalid = || 0.0;
    let value = match code {
        2 => beta_pq((x + 1.0) / 2.0, a / 2.0, a / 2.0).0,
        3 => student_pq(x, a).0,
        4 => {
            if x <= 0.0 {
                0.0
            } else {
                beta_pq(a * x / (b + a * x), a / 2.0, b / 2.0).0
            }
        }
        5 => normal_q(-x),
        6 => {
            if x <= 0.0 {
                0.0
            } else {
                gamma_pq(a / 2.0, x / 2.0).0
            }
        }
        7 => beta_pq(x, a, b).0,
        8 => {
            if x < 0.0 || a <= 0.0 || !(0.0..=1.0).contains(&b) {
                0.0
            } else if x >= a {
                1.0
            } else {
                beta_pq(1.0 - b, a - x, x + 1.0).0
            }
        }
        9 => {
            if x <= 0.0 {
                0.0
            } else {
                gamma_pq(a, x * b).0
            }
        }
        10 => {
            if x < 0.0 {
                0.0
            } else {
                gamma_pq(x + 1.0, a).1
            }
        }
        11 => {
            if b > 0.0 {
                normal_q(-(x - a) / b)
            } else {
                invalid()
            }
        }
        12 => {
            if x <= 0.0 {
                0.0
            } else {
                mixture(c / 2.0, |j| {
                    beta_pq(a * x / (b + a * x), a / 2.0 + j as f64, b / 2.0).0
                })
            }
        }
        13 => {
            if x <= 0.0 {
                0.0
            } else {
                mixture(b / 2.0, |j| gamma_pq(a / 2.0 + j as f64, x / 2.0).0)
            }
        }
        14 => {
            if b > 0.0 {
                1.0 / (1.0 + (-(x - a) / b).exp())
            } else {
                invalid()
            }
        }
        15 => {
            if b > 0.0 {
                let z = (x - a) / b;
                if z < 0.0 {
                    0.5 * z.exp()
                } else {
                    1.0 - 0.5 * (-z).exp()
                }
            } else {
                invalid()
            }
        }
        16 => {
            if b > a {
                ((x - a) / (b - a)).clamp(0.0, 1.0)
            } else {
                invalid()
            }
        }
        17 => noncentral_t(x, a, b),
        18 => {
            if b > 0.0 && c > 0.0 {
                if x <= a {
                    0.0
                } else {
                    -(-((x - a) / b).powf(c)).exp_m1()
                }
            } else {
                invalid()
            }
        }
        19 => {
            if x <= 0.0 {
                0.0
            } else {
                gamma_pq(a / 2.0, x * x / 2.0).0
            }
        }
        20 => {
            if x > 0.0 && a > 0.0 && b > 0.0 {
                let z = (b / x).sqrt();
                normal_q(-z * (x / a - 1.0)) + (2.0 * b / a).exp() * normal_q(z * (x / a + 1.0))
            } else {
                invalid()
            }
        }
        21 => {
            if b > 0.0 {
                (-(-(x - a) / b).exp()).exp()
            } else {
                invalid()
            }
        }
        22 => {
            if (0.0..=1.0).contains(&x) {
                1.0 - x
            } else {
                0.0
            }
        }
        23 => 1.0 - (-x.abs()).exp(),
        24 => 1.0 - 10.0_f64.powf(-x.abs()),
        _ => 0.0,
    };
    // f64::clamp preserves NaN: the caller will report nonconvergence.
    value.clamp(0.0, 1.0)
}
/// Invert the requested NIfTI distribution. Finite support endpoints and
/// reference BIG sentinels are handled before the monotone root search.
fn quantile(p: f64, code: i32, a: f64, b: f64, c: f64) -> f64 {
    if !(0.0..=1.0).contains(&p) || !valid_distribution(code, a, b, c) {
        return BIG;
    }
    if code == 22 {
        return 1.0 - p;
    }
    if code == 23 {
        return if p == 1.0 { BIG } else { -(-p).ln_1p() };
    }
    if code == 24 {
        return if p == 1.0 {
            BIG
        } else {
            -(-p).ln_1p() / 10.0_f64.ln()
        };
    }
    if code == 5 {
        return if p == 0.0 {
            -BIG
        } else if p == 1.0 {
            BIG
        } else {
            -normal_q_inverse(p)
        };
    }
    let (lo, hi) = match code {
        2 => (-1.0, 1.0),
        7 => (0.0, 1.0),
        8 => (0.0, a),
        16 => (a, b),
        3 | 5 | 11 | 14 | 15 | 17 | 21 => (-1.0, 1.0),
        18 => (a, a + b.max(1.0)),
        _ => (0.0, 1.0),
    };
    if p == 0.0 {
        return if matches!(code, 3 | 11 | 14 | 15 | 17 | 21) {
            -BIG
        } else {
            lo
        };
    }
    if p == 1.0 {
        return if matches!(code, 2 | 7 | 8 | 16) {
            hi
        } else {
            BIG
        };
    }
    if !(2..=21).contains(&code) {
        return BIG;
    }
    inverse(p, lo, hi, |x| cdf(x, code, a, b, c))
}

/// AFNI FI*_T2P/P2T/T2Z wrappers round arguments/results through float, and use
/// two-sided tails for t/z/correlation. The one-sided families still convert
/// to z via half their upper-tail probability, as mri_stats.c specifies.
pub(crate) fn statistical(name: &str, args: &[f64]) -> Option<f64> {
    if name == "STAT2CDF" {
        return Some(cdf(args[0], args[1] as i32, args[2], args[3], args[4]));
    }
    if name == "CDF2STAT" {
        return Some(quantile(args[0], args[1] as i32, args[2], args[3], args[4]));
    }
    let (family, operation) = name.split_once('_')?;
    let v: Vec<_> = args.iter().map(|&x| x as f32 as f64).collect();
    let x = v[0];
    let a = v.get(1).copied().unwrap_or(0.0);
    let b = v.get(2).copied().unwrap_or(0.0);
    let c = v.get(3).copied().unwrap_or(0.0);
    let code = match family {
        "FICO" => 2,
        "FITT" => 3,
        "FIFT" => 4,
        "FIZT" => 5,
        "FICT" => 6,
        "FIBT" => 7,
        "FIBN" => 8,
        "FIGT" => 9,
        "FIPT" => 10,
        _ => return None,
    };
    // FI correlation uses (sample count, fitted count, nuisance count),
    // unlike NIfTI code 2's single degrees-of-freedom parameter. Check its
    // wrapper contract separately before calling shared probability routines.
    let valid = if code == 2 {
        a > b + c
            && if operation == "P2T" {
                b >= 1.0 && c >= 1.0
            } else {
                b + c >= 1.0
            }
    } else if code == 3 {
        a >= 1.0
    } else {
        valid_distribution(code, a, b, c)
    };
    let tail = |t: f64| -> f64 {
        if !valid {
            return 1.0;
        }
        match code {
            2 => {
                if t <= 0.0 || a <= b + c || b + c < 1.0 {
                    1.0
                } else if t >= 0.9999999 {
                    0.0
                } else {
                    beta_pq(1.0 - t * t, (a - b - c) / 2.0, b / 2.0).0
                }
            }
            3 => {
                if t <= 0.0 || a < 1.0 {
                    1.0
                } else {
                    2.0 * student_pq(t, a).1
                }
            }
            5 => 2.0 * normal_q(t.abs()),
            // cdflib reports an invalid statistic above ntrial; AFNI's FI
            // wrapper converts that status to a tail of one (not zero).
            8 if t < 0.0 || t > a => 1.0,
            7 if t > 1.0 => 1.0,
            6 => {
                if t <= 0.0 {
                    1.0
                } else {
                    gamma_pq(a / 2.0, t / 2.0).1
                }
            }
            9 => {
                if t <= 0.0 {
                    1.0
                } else {
                    gamma_pq(a, t * b).1
                }
            }
            _ => 1.0 - cdf(t, code, a, b, c),
        }
    };
    let result = match operation {
        "T2P" => {
            if x == 0.0 {
                1.0
            } else {
                tail(if matches!(code, 2 | 3 | 5) {
                    x.abs()
                } else {
                    x
                })
            }
        }
        "P2T" => {
            if x >= 0.999999 {
                0.0
            } else if x <= 0.0 {
                match code {
                    2 => 0.999,
                    4 => 999.99,
                    5 => 9.99,
                    6 | 9 | 10 => 999.9,
                    7 => 0.9999,
                    _ => 99.99,
                }
            } else if !valid {
                0.0
            } else if code == 5 {
                normal_q_inverse(x / 2.0)
            } else {
                inverse(
                    1.0 - x,
                    0.0,
                    if matches!(code, 2 | 7) { 1.0 } else { 10.0 },
                    |t| 1.0 - tail(t),
                )
            }
        }
        "T2Z" => {
            if code == 5 {
                x
            } else {
                normal_q_inverse(0.5 * tail(if matches!(code, 2 | 3) { x.abs() } else { x }))
                    * if matches!(code, 2 | 3) && x < 0.0 {
                        -1.0
                    } else {
                        1.0
                    }
            }
        }
        _ => return None,
    };
    let result = result as f32 as f64;
    result.is_finite().then_some(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn probability_tails_and_inverse_identities() {
        // Two independent identities: the two tails sum to one; inversion
        // must recover the requested probability across center and tails.
        for a in [0.1, 0.5, 2.0, 30.0, 100.0] {
            for x in [0.001, 0.5, 1.0, 5.0, 40.0, 120.0] {
                let (p, q) = gamma_pq(a, x);
                assert!((p + q - 1.0).abs() < 1e-12);
                assert!((0.0..=1.0).contains(&p));
            }
        }
        for code in [
            2, 3, 4, 5, 6, 7, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
        ] {
            let (a, b, c) = if code == 16 {
                (1.0, 4.0, 0.0)
            } else {
                (4.0, 2.0, 3.0)
            };
            for p in [1e-6, 0.01, 0.3, 0.7, 0.99, 1.0 - 1e-6] {
                let x = quantile(p, code, a, b, c);
                assert!(
                    (cdf(x, code, a, b, c) - p).abs() < 1e-9,
                    "code={code}, p={p}, x={x}"
                );
            }
        }
    }

    #[test]
    fn domain_sentinels_and_nonconvergence_are_not_confused() {
        assert_eq!(cdf(1.0, 9, -1.0, 1.0, 0.0), 0.0);
        assert_eq!(quantile(0.5, 9, -1.0, 1.0, 0.0), BIG);
        assert!(inverse(0.5, 0.0, 1.0, |_| 0.0).is_nan());
        // A valid but deliberately over-budget noncentral calculation fails;
        // it must never masquerade as an ordinary CDF of zero.
        assert!(cdf(2.0, 13, 4.0, 1e8, 0.0).is_nan());
        for (family, parameters, endpoint) in [
            ("FICO", vec![30.0, 1.0, 2.0], 0.999),
            ("FITT", vec![10.0], 99.99),
            ("FIFT", vec![3.0, 12.0], 999.99),
            ("FIZT", vec![], 9.99),
            ("FICT", vec![4.0], 999.9),
            ("FIBT", vec![2.0, 4.0], 0.9999),
            ("FIBN", vec![12.0, 0.4], 99.99),
            ("FIGT", vec![2.0, 0.5], 999.9),
            ("FIPT", vec![3.0], 999.9),
        ] {
            let mut args = vec![0.0];
            args.extend(parameters);
            assert_eq!(
                statistical(&format!("{family}_P2T"), &args),
                Some(endpoint as f32 as f64)
            );
            args[0] = 1.0;
            assert_eq!(statistical(&format!("{family}_P2T"), &args), Some(0.0));
        }
        for name in ["FIFT_T2Z", "FIBT_T2Z", "FIBN_T2Z", "FIGT_T2Z"] {
            assert!(statistical(name, &[-2.0, 3.0, 0.5]).unwrap().abs() < 1e-6);
        }
    }
}
