// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Independent scientific invariants and public-library validation.
//!
//! Expected polynomial values below are hand evaluated, not obtained by calling
//! the implementation's recurrence. Golden AFNI comparisons live separately in
//! tests/run.py baselines. unwrap/expect here mean the test itself must fail if a
//! declared valid fixture cannot be constructed; production input never uses it.

use afni_xmat_maker::{
    design::{Nuisance, NuisanceKind, build_design},
    model::{Acquisition, ColumnIndex, FullRow, RetainedRow},
    table::NumericTable,
    xmat::{RowSelection, write_matrix},
};

fn nuisance(label: &str, kind: NuisanceKind, nrows: usize, values: Vec<f64>) -> Nuisance {
    let table = NumericTable::new(nrows, values.len() / nrows, values).unwrap();
    Nuisance::new(label, kind, table).unwrap()
}

fn close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (a, e) in actual.iter().zip(expected) {
        assert!((a - e).abs() < 1e-14, "{a} != {e}");
    }
}

#[test]
fn censoring_does_not_change_trends_or_nuisance_means() {
    let acquisition = Acquisition::new(vec![5], 2.0, 2, None).unwrap();
    let source = nuisance(
        "motion",
        NuisanceKind::StimFile,
        5,
        vec![10., 11., 12., 13., 14.],
    );
    let design = build_design(
        acquisition,
        vec![source],
        Some(&[false, true, false, true, true]),
    )
    .unwrap();
    close(
        design.column(ColumnIndex(0)).unwrap(),
        &[1., 1., 1., 1., 1.],
    );
    close(
        design.column(ColumnIndex(1)).unwrap(),
        &[-1., -0.5, 0., 0.5, 1.],
    );
    // P2=(3*x*x-1)/2, with its discrete full-run mean 1/4 subtracted.
    close(
        design.column(ColumnIndex(2)).unwrap(),
        &[0.75, -0.375, -0.75, -0.375, 0.75],
    );
    close(
        design.column(ColumnIndex(3)).unwrap(),
        &[-2., -1., 0., 1., 2.],
    );
    assert_eq!(
        design.retained_rows(),
        &[FullRow(1), FullRow(3), FullRow(4)]
    );
    assert_eq!(design.full_row(RetainedRow(1)), Some(FullRow(3)));
    assert_eq!(design.full_row(RetainedRow(3)), None);
}

#[test]
fn run_local_baselines_and_global_nuisance_remain_distinct() {
    let acquisition = Acquisition::new(vec![2, 3], 0.8, 0, None).unwrap();
    let source = nuisance(
        "motion",
        NuisanceKind::StimFile,
        5,
        vec![0., 0., 5., 5., 5.],
    );
    let design = build_design(acquisition, vec![source], None).unwrap();
    close(
        design.column(ColumnIndex(0)).unwrap(),
        &[1., 1., 0., 0., 0.],
    );
    close(
        design.column(ColumnIndex(1)).unwrap(),
        &[0., 0., 1., 1., 1.],
    );
    close(
        design.column(ColumnIndex(2)).unwrap(),
        &[-3., -3., 2., 2., 2.],
    );
    assert_eq!(design.acquisition().run_starts(), &[FullRow(0), FullRow(2)]);
}

#[test]
fn wide_and_zero_designs_are_not_rejected_by_an_ols_rule() {
    let acquisition = Acquisition::new(vec![2], 1., 3, None).unwrap();
    let zero = nuisance("zero", NuisanceKind::Ortvec, 2, vec![0., 0.]);
    let design = build_design(acquisition, vec![zero], None).unwrap();
    assert_eq!(design.columns().len(), 5);
    assert_eq!(design.column(ColumnIndex(4)).unwrap(), &[0., 0.]);
}

#[test]
fn one_sample_constant_and_no_polort_defaults_are_defined() {
    let design = build_design(
        Acquisition::new(vec![1], 1., 0, None).unwrap(),
        vec![],
        None,
    )
    .unwrap();
    assert_eq!(design.value(FullRow(0), ColumnIndex(0)), Some(1.));
    assert!(Acquisition::new(vec![1], 1., 1, None).is_err());
    let source = nuisance("motion", NuisanceKind::StimFile, 2, vec![3., 5.]);
    let design = build_design(
        Acquisition::new(vec![2], 1., -1, None).unwrap(),
        vec![source],
        None,
    )
    .unwrap();
    assert_eq!(design.column(ColumnIndex(0)).unwrap(), &[3., 5.]);
}

#[test]
fn source_order_and_expanded_label_collisions_are_checked() {
    let acquisition = Acquisition::new(vec![2], 1., -1, None).unwrap();
    let a = nuisance("ort", NuisanceKind::Ortvec, 2, vec![1., 2.]);
    let b = nuisance("stim", NuisanceKind::StimFile, 2, vec![3., 4.]);
    let design = build_design(acquisition, vec![a, b], None).unwrap();
    assert_eq!(design.columns()[0].label(), "stim#0");
    assert_eq!(design.columns()[1].label(), "ort[0]#0");
    let a = nuisance("motion", NuisanceKind::Ortvec, 2, vec![1., 2.]);
    let b = nuisance("motion[0]", NuisanceKind::StimFile, 2, vec![3., 4.]);
    assert!(
        build_design(
            Acquisition::new(vec![2], 1., -1, None).unwrap(),
            vec![a, b],
            None
        )
        .is_err()
    );
}

#[test]
fn invalid_dimensions_and_censor_maps_fail_before_indexing() {
    assert!(Acquisition::new(vec![usize::MAX], 1., 0, None).is_err());
    assert!(Acquisition::new(vec![0], 1., 0, None).is_err());
    assert!(Acquisition::new(vec![2], f64::NAN, 0, None).is_err());
    assert!(NumericTable::new(2, 2, vec![1., 2.]).is_err());
    assert!(NumericTable::new(1, 1, vec![f64::INFINITY]).is_err());
    let enormous = Acquisition::new(vec![i32::MAX as usize], 1., i32::MAX - 1, None).unwrap();
    assert!(build_design(enormous, vec![], None).is_err());
    assert!(
        build_design(
            Acquisition::new(vec![3], 1., 0, None).unwrap(),
            vec![],
            Some(&[true])
        )
        .is_err()
    );
    assert!(
        build_design(
            Acquisition::new(vec![3], 1., 0, None).unwrap(),
            vec![],
            Some(&[false; 3])
        )
        .is_err()
    );
}

#[test]
fn finite_inputs_that_overflow_after_centering_are_rejected() {
    let source = nuisance(
        "extreme",
        NuisanceKind::StimFile,
        3,
        vec![-f64::MAX, -f64::MAX, f64::MAX],
    );
    assert!(
        build_design(
            Acquisition::new(vec![3], 1., 0, None).unwrap(),
            vec![source],
            None
        )
        .is_err()
    );
}

#[test]
fn writer_escapes_labels_and_preserves_full_precision() {
    let source = nuisance(
        "a&\"<b",
        NuisanceKind::StimFile,
        2,
        vec![1.2345678901234567, -1e-100],
    );
    let design = build_design(
        Acquisition::new(vec![2], 1., -1, None).unwrap(),
        vec![source],
        None,
    )
    .unwrap();
    let mut output = Vec::new();
    write_matrix(
        &design,
        RowSelection::Full,
        "command 'quoted'\nnext",
        &mut output,
    )
    .unwrap();
    let text = String::from_utf8(output).unwrap();
    assert!(text.contains("a&amp;&quot;&lt;b#0"));
    let numeric: Vec<f64> = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.parse().unwrap())
        .collect();
    assert_eq!(numeric, [1.2345678901234567, -1e-100]);
    assert!(!text.contains("\nnext"));
}
