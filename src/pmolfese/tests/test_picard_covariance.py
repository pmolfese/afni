#!/usr/bin/env python3
"""NG1b/c CLI checks against independent dense NumPy oracle and AFNI REML.

Arguments: rebuilt 3dPicard, rebuilt 3dREMLfit. Requires NumPy and nibabel.
The REML check uses its exported fitted parameters, avoiding ABfile grid
rounding assumptions. No production libraries are imported into the oracle.
"""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import nibabel as nib
import numpy as np
from numpy.testing import assert_allclose
from picard_covariance_reference import correlation, ridge_gls

PICARD, REML = map(lambda s: str(Path(s).resolve()), sys.argv[1:3])
JOBS = int(os.environ.get('PICARD_TEST_JOBS', '2'))
ENV = dict(os.environ, AFNI_FLOATSCAN='NO', AFNI_NIFTI_TYPE_WARN='NO',
           AFNI_WRITE_1D_AS_PREFIX='YES', OMP_NUM_THREADS='2')


def run(args, ok=True):
    p = subprocess.run(list(map(str, args)), env=ENV, capture_output=True, text=True, cwd=folder)
    if (p.returncode == 0) != ok:
        raise AssertionError(f'{args}\n{p.stdout}\n{p.stderr}')
    return p.stdout+p.stderr


def volume(path, y, affine=None):
    nib.save(nib.Nifti1Image(np.asarray(y.T[:, None, None, :], np.float32),
                           np.eye(4) if affine is None else affine), path)


def read(prefix):
    image = nib.load(str(prefix)+'+orig.HEAD')
    return np.asarray(image.dataobj).reshape(image.shape[0], -1).T, image.header.info


def matrix_file(path, matrix, rows, nfull, starts):
    with open(path, 'w') as f:
        f.write(f'# <matrix\n# ni_type = "3*float"\n# ni_dimen = "{len(rows)}"\n')
        for k, v in dict(ColumnLabels='base;x;y', ColumnGroups='-1,1,2',
                         GoodList=','.join(map(str, rows)), NRowFull=str(nfull),
                         RunStart=','.join(map(str, starts)), RowTR='1.5').items():
            f.write(f'# {k} = "{v}"\n')
        f.write('# >\n')
        np.savetxt(f, matrix, fmt='%.9g')
        f.write('# </matrix>\n')


with tempfile.TemporaryDirectory(prefix='picard_covariance_') as temp:
    folder = Path(temp)
    rng = np.random.default_rng(821)
    nfull = 32
    rows = np.array([i for i in range(nfull) if i not in (0, 4, 5, 19, 31)])
    runs = (rows >= 16).astype(int)
    matrix = np.c_[np.ones(len(rows)), rng.normal(size=(len(rows), 2))].astype(np.float32).astype(float)
    y = rng.normal(size=(nfull, 5)).astype(np.float32).astype(float)
    xmat, data = folder/'x.xmat.1D', folder/'y.nii'
    matrix_file(xmat, matrix, rows, nfull, [0, 16])
    volume(data, y)
    base = [PICARD, '-matrix', xmat, '-input', data, '-save_nuisance']
    cases = [('white', 0., 0.), ('ar1', .6, 0.), ('arma11', .6, .2),
             ('arma11', .8, -.5), ('arma11', -.3, .4)]
    count = 0
    for model, ar, ma in cases:
        for gap in ['preserve', 'reset']:
            opts = ['-covariance', model, '-covariance_ar', ar,
                    '-covariance_ma', ma, '-covariance_gaps', gap]
            # Fraction 1 validates GLS; rho=0 elastic path supplies a fixed
            # normalized ridge lambda, allowing a fully independent oracle.
            for penalized in [False, True]:
                out = folder/f'fit{count}'
                count += 1
                fit = ['-penalty', 'enet', '-l1_ratio', 0, '-lambda', .03,
                       '-solver_tol', 1e-10] if penalized else ['-fraction', 1]
                run(base+opts+fit+['-prefix', out, '-save_residuals', '-jobs', JOBS, '-batch_size', 2])
                betas, info = read(out)
                r = correlation(rows, runs, model, ar, ma, gap)
                for v in range(5):
                    expected = ridge_gls(matrix[:, 1:], matrix[:, :1], y[rows, v], r,
                                         len(rows)*.03 if penalized else 0)
                    assert_allclose(betas[:3, v], np.r_[expected['beta'], expected['nuisance']],
                                    atol=3e-6, rtol=3e-6)
                assert 'unit_marginal' in info['PICARD_COVARIANCE']
                assert_allclose(betas[-2:], np.tile([[ar], [ma]], (1, 5)), atol=1e-7)
    # Independently refit each original-time block fold. This detects using
    # a submatrix of the full whitening operator instead of factoring the
    # training covariance, or dropping covariance in plan_like().
    from scipy.optimize import brentq

    def fractional(x, z, response, r, fraction):
        ref = ridge_gls(x, z, response, r)
        target = fraction*np.linalg.norm(ref['beta'])
        if fraction == 1 or target == 0:
            return ref
        def delta(alpha):
            return np.linalg.norm(ridge_gls(x, z, response, r, alpha)['beta'])-target
        upper = 1.
        while delta(upper) > 0:
            upper *= 2
        return ridge_gls(x, z, response, r, brentq(delta, 0, upper, xtol=1e-13))

    cvout = folder/'predict'
    cvopts = ['-covariance', 'arma11', '-covariance_ar', .7,
              '-covariance_ma', -.3, '-fractions', '1,.8,.5',
              '-cv_predict', '-cv_block_size', 8, '-save_cv']
    run(base+cvopts+['-prefix', cvout])
    cb, ci = read(cvout)
    labels = ci['BRICK_LABS'].strip('~').split('~')
    fractions = [1., .8, .5]
    for v in range(5):
        losses = []
        for fraction in fractions:
            sse = 0.
            for fold in range(4):
                train, test = rows//8 != fold, rows//8 == fold
                r = correlation(rows[train], runs[train], 'arma11', .7, -.3)
                fit = fractional(matrix[train, 1:], matrix[train, :1],
                                 y[rows[train], v], r, fraction)
                pred = matrix[test, 1:]@fit['beta']+matrix[test, :1]@fit['nuisance']
                sse += np.sum((y[rows[test], v]-pred)**2)
            losses.append(sse/len(rows))
        winner = int(np.argmin(losses))
        assert_allclose(cb[labels.index('Picard_PredMSE'), v], losses[winner],
                        rtol=3e-6, atol=3e-6)
        assert cb[labels.index('Picard_PredIndex'), v] == winner
        fit = fractional(matrix[:, 1:], matrix[:, :1], y[rows, v],
                         correlation(rows, runs, 'arma11', .7, -.3), fractions[winner])
        assert_allclose(cb[:3, v], np.r_[fit['beta'], fit['nuisance']], atol=3e-6)
    # Imported maps, exact reuse, masks and worker/batch invariance.
    ar, ma = np.array([.6, .2, .6, .2, .8], np.float32), np.array([-.2, .1, -.2, .1, -.4], np.float32)
    amap = folder/'map.nii'
    volume(amap, np.array([ar, ma]))
    mapped = base+['-covariance', 'arma11', '-covariance_map', amap, '-fraction', .7]
    first = folder/'mapfit'
    log = run(mapped+['-prefix', first, '-batch_size', 1])
    assert 'groups=3' in log
    betas, info = read(first)
    second = folder/'mapfit2'
    run(mapped+['-prefix', second, '-batch_size', 4, '-jobs', JOBS])
    assert_allclose(read(second)[0], betas, atol=0, rtol=0)
    for v in range(5):
        out = folder/f'scalar{v}'
        run(base+['-covariance', 'arma11', '-covariance_ar', ar[v],
                  '-covariance_ma', ma[v], '-fraction', .7, '-prefix', out])
        assert_allclose(read(out)[0][:, v], betas[:, v], atol=2e-6, rtol=2e-6)
    # Validate mode must inspect the selected map, including its factorization.
    run(mapped+['-validate_only'])
    bad = folder/'bad.nii'
    volume(bad, np.array([ar, [0, 0, 1.2, 0, 0]]))
    run(base+['-covariance', 'arma11', '-covariance_map', bad, '-fraction', 1,
              '-prefix', folder/'badfit'], ok=False)
    assert not (folder/'badfit+orig.HEAD').exists()
    # Native payload retains NaN when AFNI_FLOATSCAN=NO, unlike some NIfTI
    # loader conversions. Check rejection at the loaded-value boundary.
    import gzip
    native = folder/'native_map'
    run([str(Path(REML).with_name('3dcopy')), amap, native])
    payload = folder/'native_map+orig.BRIK'
    compressed = not payload.exists()
    if compressed:
        payload = Path(str(payload)+'.gz')
    raw = gzip.decompress(payload.read_bytes()) if compressed else payload.read_bytes()
    raw_values = np.frombuffer(raw, dtype=np.float32).copy()
    raw_values[2] = np.nan
    payload.write_bytes(gzip.compress(raw_values.tobytes()) if compressed else raw_values.tobytes())
    run(base+['-covariance', 'arma11', '-covariance_map', str(native)+'+orig.HEAD',
              '-fraction', 1, '-validate_only'], ok=False)
    for opts in [['-covariance', 'ar1', '-covariance_ar', 1],
                 ['-covariance', 'white', '-prewhite_ar1', 0],
                 ['-covariance', 'arma11', '-memory_mb', .001]]:
        run(base+opts+['-fraction', 1, '-prefix', folder/'invalid'], ok=False)
    wrong_affine = np.eye(4); wrong_affine[0, 3] = 5
    volume(bad, np.array([ar, ma]), wrong_affine)
    run(base+['-covariance', 'arma11', '-covariance_map', bad, '-fraction', 1,
              '-validate_only'], ok=False)
    mask = folder/'mask.nii'
    volume(mask, np.array([[1, 0, 1, 1, 1]]))
    masked = folder/'masked'
    run(mapped+['-mask', mask, '-prefix', masked])
    assert_allclose(read(masked)[0][-2:, 1], [-1, -1])
    # Actual REML estimation -> Rvar -> Picard exact GLS, not guessed grid cells.
    rv, rb = folder/'rvar', folder/'rbeta'
    run([REML, '-input', data, '-matrix', xmat, '-Grid', 3, '-CORcut', 0,
         '-Rvar', rv, '-Rbeta', rb, '-verb'])
    params, _ = read(rv)
    exact = folder/'remlimport'
    run(base+['-covariance', 'arma11', '-covariance_map', str(rv)+'+orig.HEAD',
              '-covariance_source', 'REML_full_data_reference', '-fraction', 1,
              '-prefix', exact])
    actual, _ = read(exact)
    reference, _ = read(rb)
    # Rbeta includes all source columns; Picard exports penalized then nuisance.
    assert_allclose(actual[:3], reference[[1, 2, 0]], rtol=3e-5, atol=3e-5)
    assert_allclose(actual[-2:], params[:2], atol=0, rtol=0)
    print(f'PASS {count} scalar GLS/ridge cases; independent block CV; exact maps, failures, masks, threads; REML import parity')
