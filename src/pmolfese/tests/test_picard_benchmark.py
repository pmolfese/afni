"""Validate simulation generation and adaptive stopping without a full study."""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
import benchmark_picard_covariance as benchmark


class BenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.config = json.loads(Path(__file__).with_name(
            'picard_ng1_benchmark.json').read_text())

    def test_design_and_censoring(self):
        for smoke in [True, False]:
            g = benchmark.generator(self.config, smoke)
            other = benchmark.generator(self.config, smoke)
            assert_array_equal(g['x'], other['x'])
            assert_array_equal(g['beta'], other['beta'])
            nc = 4 if smoke else self.config['conditions']
            for r in range(g['nr']):
                trials = g['trialrun'] == r
                assert_array_equal(np.bincount(g['condition'][trials]), [2]*nc)
                outside = np.arange(g['nr']*g['nt'])//g['nt'] != r
                self.assertFalse(np.any(g['x'][outside][:, trials]))
                for condition in range(nc):
                    assert_allclose(g['compact'][:, r*nc+condition],
                                    g['x'][:, trials & (g['condition'] == condition)].sum(axis=1),
                                    atol=2e-7)
                for name, removed in [('random', round(.1*g['nt'])),
                                      ('blocks', 8 if smoke else 24)]:
                    rows = g['masks'][name]
                    self.assertTrue(np.all(np.diff(rows) > 0))
                    self.assertEqual(sum(rows//g['nt'] == r), g['nt']-removed)

    def test_stationary_noise_and_replicate_prefix(self):
        g = benchmark.generator(self.config, True)
        truth = {'model': 'ar1', 'ar': .6}
        small = benchmark.draw_noise(g, truth, 8)
        large = benchmark.draw_noise(g, truth, 10000)
        assert_allclose(small, large[:, :8], atol=2e-15)
        # Independent empirical check of stationary marginal variance, actual
        # lag-two covariance, and independence across run boundaries.
        assert_allclose(np.mean(large*large, axis=1), 1, atol=.06)
        self.assertAlmostEqual(np.mean(large[0]*large[2]), .6**2, delta=.04)
        self.assertAlmostEqual(np.mean(large[63]*large[64]), 0, delta=.04)

    def test_adaptive_stopping_and_cap(self):
        for stop_at, start in [(400, 200), (None, 1900)]:
            with tempfile.TemporaryDirectory() as work:
                output = Path(work)/'report.json'
                counts = []

                def child(command, **kwargs):
                    count = int(command[command.index('--replicates')+1])
                    counts.append(count)
                    Path(command[-1]).write_text(json.dumps({
                        'precision_target_met': count == stop_at, 'failures': []}))
                    return SimpleNamespace(returncode=0)

                args = ['benchmark', '--picard', '/unused/picard', '--remlfit',
                        '/unused/reml', '--output', str(output), '--full',
                        '--adaptive', '--censor', 'all', '--replicates', str(start)]
                with patch('sys.argv', args), patch.object(benchmark.subprocess, 'run', child):
                    benchmark.main()
                self.assertEqual(counts, [200, 400] if stop_at else [1900, 2000])
                report = json.loads(output.read_text())
                self.assertEqual(len(report['adaptive_history']), 2)
                self.assertEqual(report['precision_target_met'], stop_at is not None)


if __name__ == '__main__':
    unittest.main(verbosity=2)
