// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Analytic sampled-response anchors independent of GLMsingle's fixture values.
use afni_xmat_maker::sampled::Curve;
#[test]
fn linear_interpolation_has_zero_exterior_and_exact_knots() {
    let curve = Curve::new(vec![0.0, 2.0, 0.0], 1.0, 0.0).unwrap();
    assert_eq!(curve.sample(0.5), Some(1.0));
    assert_eq!(curve.sample(1.0), Some(2.0));
    assert_eq!(curve.sample(-0.1), Some(0.0));
    assert_eq!(curve.sample(2.1), Some(0.0));
    assert_eq!(curve.sample(f64::NAN), None);
}
#[test]
fn rectangular_integration_matches_triangle_area_without_double_duration() {
    let impulse = Curve::new(vec![0.0, 1.0, 0.0], 1.0, 0.0).unwrap();
    let response = impulse.convolve_box(1.0).unwrap();
    assert_eq!(response.values(), &[0.0, 0.5, 0.5, 0.0]);
    assert_eq!(response.duration(), 1.0);
    assert_eq!(
        impulse.convolve_box(0.0).unwrap().values(),
        impulse.values()
    );
    // Half-second duration integrates only the right-hand half of the rising
    // segment at t=1: integral(.5..1) t dt = 3/8.
    assert_eq!(impulse.convolve_box(0.5).unwrap().values()[1], 0.375);
}
#[test]
fn sampled_input_rejects_nonfinite_or_unusable_shapes() {
    for values in [vec![], vec![1.0], vec![0.0, f64::NAN]] {
        assert!(Curve::new(values, 1.0, 0.0).is_err());
    }
    assert!(Curve::new(vec![0.0, 1.0], 0.0, 0.0).is_err());
    assert!(
        Curve::new(vec![0.0, -1.0], 1.0, 0.0)
            .unwrap()
            .peak_normalized()
            .is_err()
    );
}

#[test]
fn tiny_spacing_and_normalization_overflow_are_errors_not_panics() {
    let tiny = Curve::new(vec![0.0, 1.0], f64::MIN_POSITIVE, 0.0).unwrap();
    assert!(tiny.convolve_box(1.0).is_err());
    let huge = Curve::new(vec![-f64::MAX, f64::MIN_POSITIVE], 1.0, 0.0).unwrap();
    assert!(huge.peak_normalized().is_err());
}
