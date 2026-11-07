// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! M3 checks with analytic expectations, separate from the frozen AFNI oracle.
//! `unwrap` is intentional in tests: an unexpected error should fail the test,
//! while production callers propagate it with `?` and contextual messages.
use afni_xmat_maker::{
    basis::Basis,
    design::{EventMode, EventStimulus, Source, build_with_events},
    model::{Acquisition, ColumnIndex},
    timing::{Timing, TimingMode},
};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}, tolerance={tolerance}"
    );
}

#[test]
fn interpolation_and_orthogonal_families_have_analytic_anchors() {
    let spline = Basis::parse("CSPLIN(-2,10,7)").unwrap();
    let zero = Basis::parse("CSPLINzero(-2,10,7)").unwrap();
    // Cardinal cubic splines equal the identity matrix at their knots.
    for knot in 0..7 {
        for c in 0..7 {
            near(
                spline.sample(c, -2.0 + 2.0 * knot as f64).unwrap(),
                f64::from(c == knot),
                1e-6,
            );
        }
    }
    for i in 0..49 {
        let t = -2.0 + i as f64 / 4.0;
        for c in 0..5 {
            assert_eq!(zero.sample(c, t), spline.sample(c + 1, t));
        }
    }
    let poly = Basis::parse("POLY(-2,10,20)").unwrap();
    for c in 0..20 {
        near(poly.sample(c, 10.0).unwrap(), 1.0, 1e-6);
        near(
            poly.sample(c, -2.0).unwrap(),
            if c % 2 == 0 { 1.0 } else { -1.0 },
            1e-6,
        );
        assert_eq!(poly.sample(c, 11.0), Some(0.0));
    }
    // Independent Legendre recurrence at an interior coordinate avoids
    // reusing the explicit Horner expressions used by production evaluation.
    let x = 0.25;
    let mut previous = 1.0;
    let mut current = x;
    for degree in 2..20 {
        let next = ((2 * degree - 1) as f64 * x * current - (degree - 1) as f64 * previous)
            / degree as f64;
        near(poly.sample(degree, 5.5).unwrap(), next, 2e-6);
        previous = current;
        current = next;
    }
    let sine = Basis::parse("SIN(0,4,3)").unwrap();
    near(sine.sample(0, 2.0).unwrap(), 1.0, 1e-6);
    near(sine.sample(1, 2.0).unwrap(), 0.0, 1e-6);
    near(sine.sample(2, 2.0).unwrap(), -1.0, 1e-6);
}

#[test]
fn aliases_signs_and_zero_duration_rules_are_explicit() {
    for (left, right) in [
        ("SPMG", "SPMG2"),
        ("IGFUN(2)", "BLOCK4(2)"),
        ("UBLOCK(2)", "UBLOCK4(2)"),
        ("dmBLOCK", "dmBLOCK4"),
        ("dmUBLOCK", "dmUBLOCK4"),
    ] {
        let a = Basis::parse(left).unwrap();
        let b = Basis::parse(right).unwrap();
        assert_eq!(a.components(), b.components());
        for i in 0..100 {
            let t = i as f64 * 0.5;
            if a.requires_duration() {
                assert_eq!(
                    a.sample_with_duration(t, 2.0),
                    b.sample_with_duration(t, 2.0)
                );
            } else {
                for c in 0..a.components() {
                    assert_eq!(a.sample(c, t), b.sample(c, t));
                }
            }
        }
    }
    for suffix in ["", "(2)"] {
        let a = Basis::parse(&format!("MION{suffix}")).unwrap();
        let b = Basis::parse(&format!("MIONN{suffix}")).unwrap();
        for i in 0..120 {
            near(
                a.sample(0, i as f64 / 2.0).unwrap(),
                -b.sample(0, i as f64 / 2.0).unwrap(),
                1e-6,
            );
        }
    }
    for name in ["dmBLOCK", "dmUBLOCK", "dmBLOCK5", "dmUBLOCK5"] {
        let b = Basis::parse(name).unwrap();
        assert_eq!(b.sample(0, 1.0), None); // missing duration cannot mean zero
        assert_eq!(b.sample_with_duration(1.0, 0.0), Some(0.0));
        assert_eq!(b.sample_with_duration(1.0, -1.0), None);
        assert_eq!(b.sample_with_duration(1.0, 1000.0), None);
        assert_eq!(b.sample_with_duration(f64::MAX, 1.0), None);
    }
    assert_eq!(
        Basis::parse("GAM").unwrap().sample_with_duration(1.0, 2.0),
        None
    );
}

#[test]
fn duration_normalization_and_peak_width_have_independent_meaning() {
    // dmUBLOCK(-4) means normalize against the four-second response, while
    // dmBLOCK(-4) leaves the natural amplitude. A negative peak is NOT a sign.
    let b = Basis::parse("dmUBLOCK(-4)").unwrap();
    let peak = 4.0_f64 / (-(-1.0_f64).exp_m1());
    near(b.sample_with_duration(peak, 4.0).unwrap(), 1.0, 2e-6);
    let natural = Basis::parse("dmBLOCK(-4)").unwrap();
    let fixed = Basis::parse("BLOCK(4)").unwrap();
    for i in 0..50 {
        assert_eq!(
            natural.sample_with_duration(i as f64 / 2.0, 4.0),
            fixed.sample(0, i as f64 / 2.0)
        );
    }
    let gamma = Basis::parse("GAMpw(5,4)").unwrap();
    near(gamma.sample(0, 5.0).unwrap(), 1.0, 1e-6);
    let mut roots = Vec::new();
    let mut last = gamma.sample(0, 0.0).unwrap() - 0.5;
    for i in 1..20000 {
        let t = i as f64 / 1000.0;
        let next = gamma.sample(0, t).unwrap() - 0.5;
        if last * next < 0.0 {
            roots.push(t);
        }
        last = next;
    }
    assert_eq!(roots.len(), 2);
    near(roots[1] - roots[0], 4.0, 0.005);
    // Duration convolution scales each derivative component independently.
    let spm = Basis::parse("SPMG3(2)").unwrap();
    for c in 0..3 {
        let max = (0..4000)
            .map(|i| spm.sample(c, i as f64 / 100.0).unwrap().abs())
            .fold(0.0_f64, f64::max);
        near(max, 1.0, 2e-5);
    }
}

#[test]
fn expression_coordinates_precedence_and_domain_recovery_are_analytic() {
    let b = Basis::parse("EXPR(-2,6) t~x~z~1~-t^2~2^3^x~sqrt(-t)~t/0~(-2)^x~a+t").unwrap();
    for (c, want) in [
        2.0,
        0.5,
        0.0,
        1.0,
        -4.0,
        2.0_f64.powf(3.0_f64.sqrt()),
        2.0_f64.sqrt(),
        0.0,
        -2.0,
        2.0,
    ]
    .into_iter()
    .enumerate()
    {
        near(b.sample(c, 2.0).unwrap(), want, 2e-6);
    }
    assert_eq!(b.sample(0, -3.0), Some(0.0));
    assert_eq!(b.sample(0, 7.0), Some(0.0));
    assert_eq!(b.sample(99, 1.0), None);
    // Reference domain conventions intentionally differ from Rust's defaults.
    let b = Basis::parse("EXPR(0,2) log(t)~asin(t+2)~exp(100+t)~0^t~extreme(-8,-3,t)").unwrap();
    assert_eq!(b.sample(0, 0.0), Some(0.0));
    assert_eq!(b.sample(1, 0.0), Some(2.0));
    assert!(b.sample(2, 0.0).unwrap().is_finite());
    assert_eq!(b.sample(3, 0.0), Some(0.0));
    assert_eq!(b.sample(4, 2.0), Some(2.0));
}

#[test]
fn random_streams_reset_and_match_distribution_moments() {
    // Fixed seeds make this a deterministic test, not a flaky random trial.
    // Broad moment bounds test the distributions, not the implementation's
    // recurrence. URAN/IRAN also have strict support/integrality requirements.
    for (expression, mean, variance) in [
        ("uran(2)+0*t", 1.0, 1.0 / 3.0),
        ("iran(4)+0*t", 2.0, 2.0),
        ("eran(2)+0*t", 2.0, 4.0),
        ("lran(2)+0*t", 0.0, 4.0 * std::f64::consts::PI.powi(2) / 3.0),
        ("gran(3,2)+0*t", 3.0, 4.0),
    ] {
        let mut b = Basis::parse(&format!("EXPR(0,1) {expression}")).unwrap();
        b.set_expression_seed(12345);
        let values: Vec<_> = (0..40000).map(|_| b.sample(0, 0.5).unwrap()).collect();
        b.set_expression_seed(12345);
        for &v in values.iter().take(100) {
            assert_eq!(b.sample(0, 0.5), Some(v));
        }
        b.set_expression_seed(12346);
        assert_ne!(b.sample(0, 0.5), Some(values[0]));
        let average = values.iter().sum::<f64>() / values.len() as f64;
        let var = values.iter().map(|v| (v - average).powi(2)).sum::<f64>() / values.len() as f64;
        near(average, mean, 0.06);
        near(var, variance, 0.25);
        if expression.starts_with("uran") {
            assert!(values.iter().all(|&v| (0.0..2.0).contains(&v)));
        }
        if expression.starts_with("iran") {
            assert!(
                values
                    .iter()
                    .all(|&v| (0.0..=4.0).contains(&v) && v.fract() == 0.0)
            );
        }
    }
    let b = Basis::parse("EXPR(0,1) uran(t)~uran(t)").unwrap();
    assert_ne!(b.sample(0, 0.5), b.sample(1, 0.5));
}

#[test]
fn new_multicolumn_families_preserve_run_im_and_censor_identity() {
    for formula in [
        "SPMG3(2)",
        "CSPLINzero(-2,8,6)",
        "TWOGAM(4,1,0.2,8,1,2)",
        "WAV(2)",
        "POLY(-2,8,5)",
        "TRIG(-2,8,5)",
        "EXPR(-2,8) t~x~z",
        "EXPR(-2,8) uran(t)",
    ] {
        let make = |censor: bool| {
            let acq = Acquisition::new(vec![16, 16], 1.0, -1, None).unwrap();
            let times = Timing::parse("12 2 2\n*\n", "test", TimingMode::Local, &acq).unwrap();
            let source = Source::Events(
                EventStimulus::new(
                    "event",
                    Basis::parse(formula).unwrap(),
                    times,
                    EventMode::Individual,
                )
                .unwrap(),
            );
            let keep: Vec<_> = (0..32).map(|r| r % 3 != 0).collect();
            build_with_events(acq, vec![source], if censor { Some(&keep) } else { None }).unwrap()
        };
        let full = make(false);
        let retained = make(true);
        for c in 0..full.columns().len() {
            let values = full.column(ColumnIndex(c)).unwrap();
            assert_eq!(Some(values), retained.column(ColumnIndex(c)));
            assert!(values[16..].iter().all(|&v| v == 0.0)); // no leaking tail into next run
        }
        assert_eq!(full.events().len(), retained.events().len());
        for (a, b) in full.events().iter().zip(retained.events()) {
            assert_eq!(a.column, b.column);
            assert_eq!(a.event, b.event);
            assert_eq!(a.run, b.run);
        }
    }
}

#[test]
fn every_family_rejects_malformed_parameters_without_fallback() {
    for formula in [
        "GAMpw(0,4)",
        "GAMpw(4,-1)",
        "GAM(4,1,-1)",
        "TWOGAM",
        "TWOGAM(4,1,2,8,1)",
        "TWOGAMpw(4,0,0.2,8,2)",
        "TWOGAM(4,1,0.2,8,1,-1)",
        "IGFUN(0)",
        "IGFUN5(2,-1)",
        "UBLOCK(0)",
        "UBLOCK5(1,-1)",
        "dmBLOCK(1,2)",
        "dmUBLOCK(-1000)",
        "CSPLIN(0,4,3)",
        "CSPLINzero(0,4,3)",
        "CSPLIN(4,0,5)",
        "POLY(0,4,21)",
        "SIN(0,4,0)",
        "TRIG(0,4,2)",
        "SIN(0,4,2.5)",
        "SPMG4",
        "SPMG1(-1)",
        "SPMG2(1,2)",
        "SPMG3(NaN)",
        "MION(-1)",
        "MIONN(1,2)",
        "WAV(-1)",
        "WAV(1,2,3,4,5,6,7)",
        "EXPR(2,0) t",
        "EXPR(0,2) 2",
        "EXPR(0,2) \"t + x\"",
        "EXPR(0,2) t+",
        "EXPR(0,2) sin()",
        "EXPR(0,2) sqrt(t,1)",
        "EXPR(0,2) within(t,1)",
        "EXPR(0,2) pairmax(t,1,2)",
        "EXPR(0,2) unknown(t)",
        "GAM(1e38,1e38,1)",
        "SPMG3(999)",
    ] {
        assert!(Basis::parse(formula).is_err(), "accepted {formula}");
    }
    for name in [
        "AI",
        "DAI",
        "BI",
        "DBI",
        "I0",
        "I1",
        "K0",
        "K1",
        "ABSEXTREME",
    ] {
        let error = Basis::parse(&format!("EXPR(0,2) {name}(t)")).unwrap_err();
        assert!(error.to_string().contains("unimplemented"));
    }
    let deeply_nested = format!("EXPR(0,2) {}t{}", "(".repeat(100), ")".repeat(100));
    assert!(Basis::parse(&deeply_nested).is_err());
    let a = Acquisition::new(vec![20], 1.0, 0, None).unwrap();
    let t = Timing::parse("1", "test", TimingMode::Local, &a).unwrap();
    assert!(
        EventStimulus::new(
            "event",
            Basis::parse("dmBLOCK").unwrap(),
            t,
            EventMode::Ordinary
        )
        .is_err()
    );
}
