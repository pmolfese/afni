// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Independent analytic and identity tests for the M2 event path.
//! Expected hats, onset/run identities, derivative relations, and numerical
//! quadrature are separate oracles; they do not call a second copy of Maker.
use afni_xmat_maker::{
    basis::Basis,
    design::{EventMode, EventStimulus, Source, build_with_events},
    manifest::write_manifest,
    model::{Acquisition, ColumnIndex, FullRow},
    timing::{Timing, TimingMode},
};

fn event(acq: &Acquisition, text: &str, formula: &str, mode: EventMode) -> Source {
    Source::Events(
        EventStimulus::new(
            "event",
            Basis::parse(formula).unwrap(),
            Timing::parse(text, "test", TimingMode::Local, acq).unwrap(),
            mode,
        )
        .unwrap(),
    )
}

#[test]
fn analytic_tents_and_censoring_keep_original_coordinates() {
    // At onset 2, centers 0/2/4 become acquisition times 2/4/6.
    // Removing row 3 must not shift these peaks or event metadata.
    let a = Acquisition::new(vec![10], 1.0, -1, None).unwrap();
    let source = event(&a, "2\n", "TENT(0,4,3)", EventMode::Individual);
    let mut keep = vec![true; 10];
    keep[3] = false;
    let d = build_with_events(a, vec![source], Some(&keep)).unwrap();
    for row in 0..10 {
        for (c, center) in [2.0_f64, 4.0, 6.0].iter().enumerate() {
            let expected = if (2..=6).contains(&row) {
                (1.0 - (row as f64 - center).abs() / 2.0).max(0.0)
            } else {
                0.0
            };
            assert_eq!(d.value(FullRow(row), ColumnIndex(c)).unwrap(), expected);
        }
    }
    assert_eq!(d.events().len(), 3);
    assert_eq!(d.events()[2].column.0, 2);
    assert_eq!(d.events()[2].onset, 2.0);
    assert_eq!(d.retained_rows()[3].0, 4);
}

#[test]
fn im_duplicates_are_distinct_and_ordinary_responses_sum() {
    let a = Acquisition::new(vec![30], 1.0, 0, None).unwrap();
    let source = event(&a, "12 2 2\n", "GAM", EventMode::Individual);
    let d = build_with_events(a, vec![source], None).unwrap();
    assert_eq!(
        d.events().iter().map(|m| m.onset).collect::<Vec<_>>(),
        vec![12.0, 2.0, 2.0]
    );
    assert_eq!(d.column(ColumnIndex(2)), d.column(ColumnIndex(3)));
    let a = Acquisition::new(vec![30], 1.0, 0, None).unwrap();
    let source = event(&a, "12 2 2\n", "GAM", EventMode::Ordinary);
    let ordinary = build_with_events(a, vec![source], None).unwrap();
    for r in 0..30 {
        let expected = (1..4).fold(0.0_f32, |sum, c| {
            sum + d.value(FullRow(r), ColumnIndex(c)).unwrap() as f32
        });
        assert_eq!(
            ordinary.value(FullRow(r), ColumnIndex(1)).unwrap(),
            expected as f64
        );
    }
    assert!(ordinary.events().iter().all(|m| m.column.0 == 1));
}

#[test]
fn global_column_order_and_stable_event_order_are_distinct() {
    // The second run's token occurs first globally. AFNI IM column order
    // follows tokens; stable IDs still follow run then source order.
    let a = Acquisition::new(vec![20, 20], 1.0, 0, None).unwrap();
    let timing = Timing::parse("23 4 2\n", "global", TimingMode::Global, &a).unwrap();
    let source = Source::Events(
        EventStimulus::new(
            "trial",
            Basis::parse("SPMG2").unwrap(),
            timing,
            EventMode::Individual,
        )
        .unwrap(),
    );
    let d = build_with_events(a, vec![source], None).unwrap();
    let relations: Vec<_> = d
        .events()
        .iter()
        .map(|m| (m.event.0, m.run.0, m.onset, m.basis.0, m.column.0))
        .collect();
    assert_eq!(
        relations,
        vec![
            (0, 0, 4.0, 0, 4),
            (0, 0, 4.0, 1, 5),
            (1, 0, 2.0, 0, 6),
            (1, 0, 2.0, 1, 7),
            (2, 1, 3.0, 0, 2),
            (2, 1, 3.0, 1, 3)
        ]
    );
    for row in 0..20 {
        assert_eq!(d.value(FullRow(row), ColumnIndex(2)), Some(0.0));
    }
    for row in 20..40 {
        assert_eq!(d.value(FullRow(row), ColumnIndex(4)), Some(0.0));
    }
}

#[test]
fn empty_and_censored_runs_preserve_trial_identity() {
    let a = Acquisition::new(vec![20, 20], 1.0, 0, None).unwrap();
    let source = event(&a, "19 2\n*\n", "GAM", EventMode::Individual);
    let keep: Vec<_> = (0..40).map(|r| r >= 20).collect();
    let d = build_with_events(a, vec![source], Some(&keep)).unwrap();
    assert_eq!(d.events().len(), 2);
    assert!(d.column(ColumnIndex(2)).unwrap().iter().all(|&x| x == 0.0));
    for row in 20..40 {
        assert_eq!(d.value(FullRow(row), ColumnIndex(3)), Some(0.0));
    }
    let mut bytes = Vec::new();
    write_manifest(&d, &mut bytes).unwrap();
    let manifest = String::from_utf8(bytes).unwrap();
    assert!(manifest.starts_with("# 1dXmatMaker events 1\nevent\tstimulus\tcondition\trun"));
    assert_eq!(manifest.lines().count(), 4); // Header pair + both events, even without retained support.
}

#[test]
fn gamma_and_spm_conventions_have_independent_anchors() {
    let gamma = Basis::parse("GAM(4,1.2)").unwrap();
    assert!((gamma.sample(0, 4.8).unwrap() - 1.0).abs() < 2e-6);
    assert_eq!(Basis::parse("GAM").unwrap().sample(0, 12.0), Some(0.0));
    assert!(
        Basis::parse("GAM(8.6,0.547)")
            .unwrap()
            .sample(0, 12.0)
            .unwrap()
            > 0.0
    );
    let spm = Basis::parse("SPMG2").unwrap();
    for t in [2.0, 5.0, 9.0, 17.0] {
        let h = 0.01;
        let derivative =
            (spm.sample(0, t + h).unwrap() - spm.sample(0, t - h).unwrap()) / (2.0 * h);
        assert!((derivative - spm.sample(1, t).unwrap()).abs() < 1e-5);
    }
    assert_eq!(spm.sample(0, 25.0), Some(0.0));
    assert_eq!(spm.sample(2, 5.0), None);
}

#[test]
fn block_matches_independent_midpoint_quadrature_and_peak() {
    // Midpoint integration of a peak-one gamma kernel does not use the
    // exponential-polynomial antiderivative in production. Historical rounded
    // constants are included as an explicit correction to the independent integral.
    // Otherwise the BLOCK4 early branch differs from the exact integral by
    // (54.5982-exp(4))*24/256, about 4.67e-6; this is reference behavior.
    for power in [4, 5] {
        let natural = Basis::parse(&format!("BLOCK{power}(2.5)")).unwrap();
        for t in [1.0_f64, 3.0, 5.0, 12.0] {
            let width = t.min(2.5);
            let steps = 10000;
            let dt = width / steps as f64;
            let expected = (0..steps)
                .map(|i| {
                    let u = t - (i as f64 + 0.5) * dt;
                    (u / power as f64).powi(power) * (power as f64 - u).exp() * dt
                })
                .sum::<f64>();
            let expected = if power == 4 && t < 2.5 {
                expected + (54.5982 - 4.0_f64.exp()) * 24.0 / 256.0
            } else if power == 5 {
                expected * 148.413 / 5.0_f64.exp()
            } else {
                expected
            };
            assert!((natural.sample(0, t).unwrap() - expected).abs() < 1e-6);
        }
        let normalized = Basis::parse(&format!("BLOCK{power}(2.5,2)")).unwrap();
        let peak = 2.5 / (1.0 - (-2.5 / power as f64).exp());
        assert!((normalized.sample(0, peak).unwrap() - 2.0).abs() < 2e-6);
        assert_eq!(normalized.duration(), 2.5);
        assert_eq!(normalized.sample(0, 17.5), Some(0.0));
    }
}

#[test]
fn malformed_timing_and_basis_do_not_become_other_models() {
    let a = Acquisition::new(vec![20, 20], 1.0, 0, None).unwrap();
    for text in [
        "2\n",
        "2\n3\n4\n",
        "*\n*\n",
        "-1\n3\n",
        "NaN\n3\n",
        "20\n3\n",
        "2*3\n3\n",
        "2:3\n3\n",
        "2 *\n3\n",
    ] {
        assert!(
            Timing::parse(text, "invalid", TimingMode::Local, &a).is_err(),
            "{text}"
        );
    }
    for model in [
        "SPMG9",
        "SPMG1(-1)",
        "GAM(8.6,.547,-2)",
        "GAM(0,1)",
        "BLOCK(0)",
        "BLOCK(2,-1)",
        "TENT(0,2,2.5)",
        "TENT(0,2,2.00000001)",
        "TENTzero(0,2,2)",
        "TENT(2,0,3)",
        "GAM(NaN,2)",
        "GAM(1,2)junk",
        "EXPR(0,2)2",
    ] {
        assert!(Basis::parse(model).is_err(), "{model}");
    }
    assert!(Timing::parse("19.9999999", "boundary", TimingMode::Global, &a).is_err());
    let too_long = Acquisition::new(vec![(1 << 24) + 1], 1.0, 0, None).unwrap();
    assert!(Timing::parse("2", "large", TimingMode::Local, &too_long).is_err());
    let source = event(&a, "2\n3\n", "GAM", EventMode::Individual);
    let different = Acquisition::new(vec![20, 20], 0.8, 0, None).unwrap();
    assert!(build_with_events(different, vec![source], None).is_err());
}
