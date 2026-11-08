"""NG1a known-answer/reference checks. Run directly with Python + NumPy."""
import json
from pathlib import Path
import unittest
import numpy as np
from numpy.testing import assert_allclose
from picard_covariance_reference import correlation, ridge_gls


class CovarianceReferenceTests(unittest.TestCase):
    def test_hand_covariance_and_gaps(self):
        assert_allclose(correlation([0, 2, 3], [0, 0, 1], 'ar1', .5),
                        [[1, .25, 0], [.25, 1, 0], [0, 0, 1]])
        assert_allclose(correlation([0, 2], [0, 0], 'ar1', .5, gap='reset'), np.eye(2))
        assert_allclose(correlation([8], [3], 'arma11', .8, -.5), [[1]])

    def test_stationary_impulse_reference(self):
        # Independent convolution of MA-infinity coefficients, not the closed
        # form in correlation(). Tail truncation is negligible for |ar|<=.8.
        for ar, ma in [(0, .4), (.6, .2), (.8, -.5), (-.4, .3), (.6, -.6)]:
            psi = np.r_[1., (ar+ma)*ar**np.arange(1000)]
            acv = np.array([psi[:len(psi)-k]@psi[k:] for k in range(10)])
            actual = correlation(np.arange(10), np.zeros(10, int), 'arma11', ar, ma)
            assert_allclose(actual[0], acv/acv[0], atol=2e-14, rtol=2e-14)
            np.linalg.cholesky(actual)

    def test_nested_limits(self):
        rows, runs = [0, 1, 4, 9], [0, 0, 0, 1]
        assert_allclose(correlation(rows, runs), np.eye(4))
        assert_allclose(correlation(rows, runs, 'arma11', .7, 0),
                        correlation(rows, runs, 'ar1', .7))
        assert_allclose(correlation(rows, runs, 'arma11', .7, -.7), np.eye(4))

    def test_fold_submatrix(self):
        full = correlation(np.arange(7), np.zeros(7, int), 'arma11', .5, -.2)
        keep = [0, 2, 6]
        assert_allclose(correlation(keep, [0]*3, 'arma11', .5, -.2),
                        full[np.ix_(keep, keep)])
        lower = np.linalg.cholesky(full)
        white = np.linalg.solve(lower, full)
        assert_allclose(np.linalg.solve(lower, white.T).T, np.eye(7), atol=1e-14)

    def test_hand_ridge(self):
        z = np.empty((2, 0))
        assert_allclose(ridge_gls(np.eye(2), z, [2, 4], np.eye(2), 2)['beta'], [2/3, 4/3])
        r = correlation([0, 1], [0, 0], 'ar1', .5)
        # X=I, penalty=1: beta=(I+R)^-1 y, an explicit 2x2 answer.
        assert_allclose(ridge_gls(np.eye(2), z, [2, 4], r, 1)['beta'], [8/15, 28/15])

    def test_nuisance_and_rank(self):
        x = np.array([[-1.], [0.], [1.]])
        z = np.ones((3, 2))
        out = ridge_gls(x, z, [1, 2, 3], np.eye(3), 2)
        assert_allclose(out['beta'], [.5])
        assert_allclose(out['nuisance'], [1, 1])
        self.assertEqual((out['xrank'], out['zrank']), (1, 1))
        redundant = ridge_gls(np.c_[x, x], z, [1, 2, 3], np.eye(3))
        assert_allclose(redundant['beta'], [.5, .5])
        explained = ridge_gls(z, z, [2, 2, 2], np.eye(3))
        assert_allclose(explained['beta'], [0, 0])

    def test_joint_normal_equation_reference(self):
        rng = np.random.default_rng(21)
        for n, p in [(8, 3), (4, 7)]:
            x, z, y = rng.normal(size=(n, p)), np.ones((n, 1)), rng.normal(size=n)
            r = correlation(np.arange(n), np.zeros(n, int), 'arma11', .6, -.3)
            d = np.c_[x, z]
            expected = np.linalg.solve(d.T@np.linalg.solve(r, d)
                                       + np.diag([.7]*p+[0]), d.T@np.linalg.solve(r, y))
            out = ridge_gls(x, z, y, r, .7)
            assert_allclose(np.r_[out['beta'], out['nuisance']], expected, atol=1e-12)
            assert_allclose(out['fitted']+out['residual'], y, atol=1e-14)

    def test_legacy_scale_conversion(self):
        rho, penalty = .6, .8
        rows = [0, 1, 4, 5]
        r = correlation(rows, [0]*4, 'ar1', rho, gap='reset')
        old = np.diag([np.sqrt(1-rho*rho), 1]*2)
        old[1, 0] = old[3, 2] = -rho
        x, y = np.array([[1.], [2.], [-1.], [3.]]), np.array([2., 1., 0., 4.])
        expected = np.linalg.solve((old@x).T@(old@x)+penalty*np.eye(1),
                                   (old@x).T@(old@y))
        out = ridge_gls(x, np.empty((4, 0)), y, r, penalty/(1-rho*rho))
        assert_allclose(out['beta'], expected, atol=1e-14)

    def test_invalid_inputs(self):
        for rows, runs in [([], []), ([1, 1], [0, 0]), ([0, 1], [1, 0]),
                           ([0., 1.], [0, 0]), ([-1], [0])]:
            with self.assertRaises(ValueError): correlation(rows, runs)
        for kw in [dict(ar=1, model='ar1'), dict(ma=np.nan), dict(gap='fill'),
                   dict(ar=.2), dict(model='ar1', ma=.1)]:
            with self.assertRaises(ValueError): correlation([0], [0], **kw)
        with self.assertRaises(ValueError):
            ridge_gls(np.eye(2), np.empty((2, 0)), [1, 2], np.eye(2), -1)
        with self.assertRaises(np.linalg.LinAlgError):
            ridge_gls(np.eye(2), np.empty((2, 0)), [1, 2], [[1, 2], [2, 1]])

    def test_benchmark_contract(self):
        cfg = json.loads(Path(__file__).with_name('picard_ng1_benchmark.json').read_text())
        self.assertEqual(cfg['fitted_models'], ['white', 'ar1', 'arma11'])
        self.assertEqual(len(cfg['outer_evaluation_runs']), 2)
        self.assertFalse(set(cfg['tuning_runs']) & set(cfg['outer_evaluation_runs']))
        self.assertEqual(cfg['status'], 'frozen_protocol_not_executed')
        self.assertTrue(all(correlation([0, 1], [0, 0], **c).shape == (2, 2)
                            for c in cfg['noise_cases']))


if __name__ == '__main__':
    unittest.main(verbosity=2)
