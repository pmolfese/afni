// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Independent M4 invariants. These tests use analytic triangular responses
//! and explicit coefficient vectors, rather than duplicating the builder.
use afni_xmat_maker::{
    basis::Basis,
    censor,
    design::{EventMode, EventStimulus, Source, build_with_events},
    input,
    model::{Acquisition, ColumnIndex, FullRow},
    timing::{Timing, TimingFormat, TimingMode, TimingOptions},
    xmat::{RowSelection, write_matrix},
};

fn am(mode: EventMode, raw: bool) -> afni_xmat_maker::design::Design {
    let a = Acquisition::new(vec![6, 6], 1.0, 0, None).unwrap();
    let t = Timing::parse_with_options(
        "0*1,3 3*5,7\n0*9,11",
        "test",
        TimingMode::Local,
        &a,
        TimingOptions {
            format: TimingFormat::Married { duration: false },
            ..Default::default()
        },
    )
    .unwrap();
    let s = EventStimulus::new("trial", Basis::parse("TENT(0,2,3)").unwrap(), t, mode)
        .unwrap()
        .centered(raw, None)
        .unwrap();
    let mut keep = vec![true; 12];
    keep[0] = false;
    keep[6] = false;
    build_with_events(a, vec![Source::Events(s)], Some(&keep)).unwrap()
}
#[test]
fn am2_centers_before_censoring_and_expands_each_component() {
    let d = am(EventMode::Am2, false);
    let raw = am(EventMode::Am2, true);
    assert_eq!(d.columns().len(), 11);
    assert_eq!(d.events().len(), 27);
    // Means are (5,7), including the two events whose onset rows are removed.
    for r in 0..12 {
        for c in 0..3 {
            let v = |design: &afni_xmat_maker::design::Design, c| {
                design.value(FullRow(r), ColumnIndex(c)).unwrap()
            };
            assert_eq!(v(&raw, 5 + c), v(&d, 5 + c) + 5.0 * v(&d, 2 + c));
            assert_eq!(v(&raw, 8 + c), v(&d, 8 + c) + 7.0 * v(&d, 2 + c));
        }
    }
    assert_eq!(d.value(FullRow(0), ColumnIndex(5)), Some(-4.0));
    assert_eq!(d.value(FullRow(6), ColumnIndex(5)), Some(4.0));
    assert_eq!(d.retained_rows().len(), 10);
}
#[test]
fn am1_sums_modulators_without_centering() {
    let d = am(EventMode::Am1, false);
    assert_eq!(d.columns().len(), 5);
    for (row, amp) in [(0, 4.0), (3, 12.0), (6, 20.0)] {
        for component in 0..3 {
            assert_eq!(
                d.value(FullRow(row + component), ColumnIndex(2 + component)),
                Some(amp)
            );
        }
    }
}
#[test]
fn contrast_rows_select_final_columns_and_repeated_terms_assign() {
    let mut d = am(EventMode::Am2, false);
    d.add_symbolic_contrast("effects", "SYM: +trial[[0..1]] -trial[[3..4]] +Ort[0]")
        .unwrap();
    d.add_symbolic_contrast("last", "SYM: 2*trial[*] -trial[1]")
        .unwrap();
    let mut out = Vec::new();
    write_matrix(&d, RowSelection::Retained, "test", &mut out).unwrap();
    let text = String::from_utf8(out).unwrap();
    assert!(
        text.contains("GltMatrix_000000 = \"2,11,1,0,1,0,0,-1,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,0,0\"")
    );
    assert!(text.contains("GltMatrix_000001 = \"1,11,0,0,2,-1,2,2,2,2,2,2,2\""));
    for invalid in [
        "SYM: unknown",
        "SYM: trial[99]",
        "SYM: trial[[0..1]] trial[[2..4]]",
        "SYM: trial +",
    ] {
        assert!(d.add_symbolic_contrast("bad", invalid).is_err());
    }
}
#[test]
fn censor_rules_keep_original_global_and_run_coordinates() {
    let a = Acquisition::new(vec![6, 8], 1.0, 0, None).unwrap();
    let keep = censor::combine(&a, vec![true; 14], 1, Some(5), &["1:2,4 *:3".into()]).unwrap();
    assert_eq!(
        keep.iter()
            .enumerate()
            .filter_map(|(i, &v)| v.then_some(i))
            .collect::<Vec<_>>(),
        vec![1, 5, 7, 8, 10, 11]
    );
    assert!(censor::combine(&a, vec![true; 14], 0, None, &["2:8".into()]).is_err());
}
#[test]
fn selector_order_steps_repetition_and_bounds_are_explicit() {
    assert_eq!(
        input::indices("$..0(2),1,1", 6).unwrap(),
        vec![5, 3, 1, 1, 1]
    );
    assert_eq!(input::indices("1..$-1", 5).unwrap(), vec![1, 2, 3]);
    for s in ["6", "0..6", "$-9", "0..2(0)", ""] {
        assert!(input::indices(s, 6).is_err());
    }
}
#[test]
fn invalid_married_fields_and_empty_im_never_become_zeros() {
    let a = Acquisition::new(vec![20], 1.0, 0, None).unwrap();
    let opts = TimingOptions {
        format: TimingFormat::Married { duration: true },
        allow_empty: true,
        ..Default::default()
    };
    for text in [
        "2*NaN:1",
        "2*1:1000",
        "2*1:-1",
        "2*1",
        "2:1 4*1:2",
        "-2*1:2",
        "20*1:2",
    ] {
        assert!(
            Timing::parse_with_options(text, "test", TimingMode::Local, &a, opts).is_err(),
            "{text}"
        );
    }
    let t = Timing::parse_with_options("*", "test", TimingMode::Local, &a, opts).unwrap();
    assert!(
        EventStimulus::new(
            "t",
            Basis::parse("dmBLOCK").unwrap(),
            t,
            EventMode::Individual
        )
        .is_err()
    );
}

#[test]
fn expanded_zero_glt_terms_fall_back_to_ordinary_coefficients() {
    let mut d = am(EventMode::Am2, false);
    d.add_symbolic_contrast("zero", "SYM: 2*trial[0] +trial[[0]] 0*trial[[0]]")
        .unwrap();
    let mut out = Vec::new();
    write_matrix(&d, RowSelection::Retained, "test", &mut out).unwrap();
    assert!(
        String::from_utf8(out)
            .unwrap()
            .contains("GltMatrix_000000 = \"1,11,0,0,2,0,0,0,0,0,0,0,0\"")
    );
}
