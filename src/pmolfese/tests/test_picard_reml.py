#!/usr/bin/env python3
"""Native REML likelihood, rank, subset isolation and failure tests."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import nibabel as nib
import numpy as np
from numpy.testing import assert_allclose
from picard_covariance_reference import correlation, reml_grid, ridge_gls

BINARY = str(Path(sys.argv[1]).resolve())
ENV = dict(os.environ, AFNI_FLOATSCAN='NO', AFNI_WRITE_1D_AS_PREFIX='YES',
           OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')


def xmat(path, matrix, rows, nfull, starts):
    n, p = matrix.shape
    with open(path, 'w') as f:
        f.write(f'# <matrix\n# ni_type = "{p}*float"\n# ni_dimen = "{n}"\n')
        for k, val in {'ColumnLabels': ';'.join(f'x{i}' for i in range(p)),
                       'ColumnGroups': ','.join(['-1']+['1']*(p-1)),
                       'GoodList': ','.join(map(str, rows)), 'NRowFull': nfull,
                       'RunStart': ','.join(map(str, starts)), 'RowTR': 1.5}.items():
            f.write(f'# {k} = "{val}"\n')
        f.write('# >\n'); np.savetxt(f, matrix, fmt='%.9g'); f.write('# </matrix>\n')


def volume(path, y):
    nib.save(nib.Nifti1Image(y.T[:, None, None, :].astype(np.float32), np.eye(4)), path)


def read(path):
    img = nib.load(str(path)+'+orig.HEAD')
    return np.asarray(img.dataobj).reshape(img.shape[0], -1).T, img.header.info


def run(args, folder, ok=True):
    p = subprocess.run([BINARY, *map(str, args)], env=ENV, cwd=folder,
                       capture_output=True, text=True)
    assert (p.returncode == 0) == ok, p.stderr
    return p


def main():
    with tempfile.TemporaryDirectory(prefix='picard_reml_') as work:
        folder = Path(work)
        rng = np.random.default_rng(122)
        rows = np.array([i for i in range(30) if i not in (3, 16, 25)])
        runs = (rows >= 15).astype(int)
        design = np.c_[np.ones(len(rows)), rng.normal(size=len(rows))].astype(np.float32).astype(float)
        pilot = np.c_[design, design[:, 1]].astype(np.float32).astype(float)
        y = rng.normal(size=(30, 3)).astype(np.float32).astype(float)
        xmat(folder/'fit.1D', design, rows, 30, [0, 15])
        xmat(folder/'pilot.1D', pilot, rows, 30, [0, 15])
        volume(folder/'data.nii', y)
        selected = np.arange(18)
        np.savetxt(folder/'rows.1D', rows[selected], fmt='%d')
        base = ['-input', folder/'data.nii', '-matrix', folder/'fit.1D',
                '-covariance_estimate', folder/'pilot.1D',
                '-covariance_rows', folder/'rows.1D', '-fraction', 1]
        for model in ['white', 'ar1', 'arma11']:
            for gap in ['preserve', 'reset']:
                out = folder/(model+gap)
                opts = base+['-covariance', model, '-covariance_gaps', gap]
                run(opts+['-prefix', out], folder)
                actual, info = read(out)
                expected, _ = reml_grid(pilot[selected], y[rows[selected]],
                                        rows[selected], runs[selected], model, gap)
                assert_allclose(actual[-6], expected[:, 2], rtol=4e-5, atol=4e-5)
                assert_allclose(actual[-5], expected[:, 3], rtol=4e-5, atol=4e-5)
                assert_allclose(actual[-4], 2)
                for v in range(3):
                    ar, ma = actual[-2:, v]
                    assert_allclose(correlation(rows, runs, model, ar, ma, gap),
                                    correlation(rows, runs, model, *expected[v, :2], gap), atol=1e-6)
                    fit = ridge_gls(design[:, 1:], design[:, :1], y[rows, v],
                                    correlation(rows, runs, model, ar, ma, gap))
                    assert_allclose(actual[0, v], fit['beta'][0], atol=3e-6)
                assert 'native_grid_REML' in info['PICARD_COVARIANCE_ESTIMATION']
        # Batch changes preserve every learned parameter and fitted value.
        run(base+['-covariance', 'arma11', '-batch_size', 1,
                  '-prefix', folder/'batch'], folder)
        assert_allclose(read(folder/'batch')[0], read(folder/'arma11preserve')[0],
                        atol=0, rtol=0)
        # Independent singleton reset segments carry no correlation information.
        np.savetxt(folder/'isolated.1D', [0, 2, 6, 8, 10, 12], fmt='%d')
        run(base+['-covariance', 'arma11', '-covariance_gaps', 'reset',
                  '-covariance_rows', folder/'isolated.1D',
                  '-prefix', folder/'flat'], folder)
        flat = read(folder/'flat')[0]
        assert_allclose(flat[-3], 3)
        assert_allclose(flat[-2:], 0)
        # Estimation cannot see held-out response changes.
        y[rows[18:]] += 100
        volume(folder/'changed.nii', y)
        changed = [folder/'changed.nii' if value == folder/'data.nii' else value for value in base]
        run(changed+['-covariance', 'arma11', '-prefix', folder/'changed'], folder)
        assert_allclose(read(folder/'changed')[0][-6:], read(folder/'arma11preserve')[0][-6:], atol=0, rtol=0)
        xmat(folder/'saturated.1D', np.eye(len(rows)), rows, 30, [0, 15])
        bad = [folder/'saturated.1D' if value == folder/'pilot.1D' else value for value in base]
        run(bad+['-covariance', 'arma11', '-prefix', folder/'bad'], folder, False)
        volume(folder/'zero.nii', np.zeros((30, 3)))
        bad = [folder/'zero.nii' if value == folder/'data.nii' else value for value in base]
        run(bad+['-covariance', 'arma11', '-prefix', folder/'bad'], folder, False)
        rejected = run(base[:-2]+['-covariance', 'ar1', '-fractions', '1,.5',
                        '-cv_predict', '-prefix', folder/'bad'], folder, False)
        assert 'native covariance requires' in rejected.stderr
        for invalid_rows in [[0, 0], [0, 3], [0, 30], [0, 1.5]]:
            np.savetxt(folder/'invalid_rows.1D', invalid_rows)
            run(base+['-covariance', 'ar1', '-covariance_rows',
                      folder/'invalid_rows.1D', '-prefix', folder/'bad'], folder, False)
        original_pilot = (folder/'pilot.1D').read_text()
        (folder/'wrong_tr.1D').write_text(original_pilot.replace(
            'RowTR = "1.5"', 'RowTR = "2.0"'))
        rejected = run(base+['-covariance', 'ar1', '-covariance_estimate',
                            folder/'wrong_tr.1D', '-prefix', folder/'bad'], folder, False)
        assert 'pilot TR mismatch' in rejected.stderr
        assert not (folder/'bad+orig.HEAD').exists()
        # Seeded parameter recovery probe: independent stationary AR1 draws,
        # not a test that every realization must identify the same grid cell.
        rr = np.arange(160)
        recovery_x = np.c_[np.ones(160), np.linspace(-1, 1, 160)]
        recovery_y = np.linalg.cholesky(correlation(rr, np.zeros(160, int), 'ar1', .6))@rng.normal(size=(160, 24))
        volume(folder/'recovery.nii', recovery_y)
        xmat(folder/'recovery.1D', recovery_x, rr, 160, [0])
        run(['-input', folder/'recovery.nii', '-matrix', folder/'recovery.1D',
             '-covariance', 'ar1', '-covariance_estimate', folder/'recovery.1D',
             '-fraction', 1, '-prefix', folder/'recovery'], folder)
        recovered = read(folder/'recovery')[0]
        assert abs(np.mean(recovered[-2])-.6) < .08
        assert_allclose(recovered[-1], 0)
        assert .8 < np.mean(recovered[-6]) < 1.2
        print('AR1 recovery: mean a=', float(np.mean(recovered[-2])),
              'mean marginal variance=', float(np.mean(recovered[-6])))
        print('PASS native white/AR1/ARMA REML, gaps, QR oracle, redundant pilot, subset isolation, failures')


if __name__ == '__main__':
    main()
