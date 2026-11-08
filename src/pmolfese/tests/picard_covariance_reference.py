"""Independent dense NG1 oracle, not a production Picard implementation.

All arrays use rows=retained time observations. Covariance uses original integer
TR indices and independent runs. NumPy/LAPACK Cholesky and least squares are
intentionally separate from AFNI's production kernels. The QR grid estimator
is a reference oracle, not a production noise-estimation implementation.
"""
import numpy as np


def correlation(rows, runs, model='white', ar=0.0, ma=0.0, gap='preserve'):
    """Return unit-diagonal stationary correlation on retained observations.

    ARMA convention: e[t] = ar*e[t-1] + u[t] + ma*u[t-1]. Parameters
    are scalar, finite and strictly inside (-1,1). Rows are increasing original
    integer TR indices; run IDs are nonnegative, nondecreasing integers.
    Reset creates independent stationary segments at every missing TR.
    Invalid/empty timelines and incompatible model parameters raise ValueError.
    """
    rows, runs = np.asarray(rows), np.asarray(runs)
    if (rows.ndim != 1 or runs.shape != rows.shape or rows.size == 0
            or rows.dtype.kind not in 'iu' or runs.dtype.kind not in 'iu'
            or np.any(rows < 0) or np.any(runs < 0)
            or np.any(rows[1:] <= rows[:-1])
            or np.any(runs[1:] < runs[:-1])):
        raise ValueError('invalid retained timeline')
    if (model not in ('white', 'ar1', 'arma11')
            or gap not in ('preserve', 'reset')
            or not np.isfinite([ar, ma]).all() or abs(ar) >= 1 or abs(ma) >= 1
            or (model == 'white' and (ar != 0 or ma != 0))
            or (model == 'ar1' and ma != 0)):
        raise ValueError('invalid covariance specification')
    # Nonnegative signed indices cannot overflow their pairwise difference.
    # Keep a Python-integer fallback for unusually large uint64 timelines.
    if int(rows.max()) <= np.iinfo(np.int64).max:
        signed = rows.astype(np.int64)
        lag = np.abs(signed[:, None]-signed[None, :])
    else:
        lag = np.array([[abs(int(i)-int(j)) for j in rows] for i in rows])
    if model == 'white':
        result = np.eye(len(rows))
    else:
        # gamma(1)/gamma(0), derived from the stationary MA-infinity expansion.
        first = (ar + ma) * (1 + ar * ma) / (1 + 2 * ar * ma + ma * ma)
        result = first * np.power(ar, np.maximum(lag-1, 0))
        np.fill_diagonal(result, 1.0)
    same = runs[:, None] == runs[None, :]
    if gap == 'reset':
        starts = np.r_[True, (rows[1:] != rows[:-1]+1)
                       | (runs[1:] != runs[:-1])]
        segments = np.cumsum(starts)
        same &= segments[:, None] == segments[None, :]
    return np.where(same, result, 0.0)


def ridge_gls(x, z, y, covariance, penalty=0.0, rank_tol=1e-10):
    """Fit ||L^-1(y-Xb-Zg)||^2 + penalty*||b||^2 in float64.

    penalty is the XtX diagonal penalty (n times normalized lambda). Z is
    unpenalized and may be redundant or have zero columns. Return minimum-norm
    projected beta, minimum-norm recovered nuisance coefficients, original-unit
    fitted/residual arrays, weighted SSE, and weighted-design ranks. y is a
    vector. All returned arrays are new. Invalid inputs raise ValueError;
    non-positive-definite covariance raises numpy.linalg.LinAlgError.
    """
    x, z, y, r = [np.asarray(a, dtype=float) for a in (x, z, y, covariance)]
    if (y.ndim != 1 or len(y) == 0 or x.ndim != 2 or z.ndim != 2
            or x.shape[0] != len(y) or z.shape[0] != len(y)
            or r.shape != (len(y), len(y))
            or not all(np.isfinite(a).all() for a in (x, z, y, r))
            or not np.isfinite(penalty) or penalty < 0
            or not np.isfinite(rank_tol) or not 0 < rank_tol < 1
            or not np.allclose(r, r.T, rtol=0, atol=1e-14)):
        raise ValueError('invalid GLS inputs')
    lower = np.linalg.cholesky(r)
    wx, wz, wy = [np.linalg.solve(lower, a) for a in (x, z, y)]
    u, s, _ = np.linalg.svd(wz, full_matrices=False)
    zrank = int(np.sum(s > rank_tol*s[0])) if s.size else 0
    q = u[:, :zrank]
    px, py = wx-q@(q.T@wx), wy-q@(q.T@wy)
    # Match the declared projection cancellation policy for fully nuisance-
    # explained columns/responses; do not normalize retained design columns.
    eps = 64*np.finfo(float).eps*max(wz.shape)
    px[:, np.linalg.norm(px, axis=0) <= eps*np.linalg.norm(wx, axis=0)] = 0
    if np.linalg.norm(py) <= eps*np.linalg.norm(wy):
        py[:] = 0
    singular = np.linalg.svd(px, compute_uv=False)
    xrank = int(np.sum(singular > rank_tol*singular[0])) if singular.size else 0
    if penalty == 0:
        beta = np.linalg.lstsq(px, py, rcond=rank_tol)[0]
    else:
        augmented = np.vstack((px, np.sqrt(penalty)*np.eye(x.shape[1])))
        response = np.r_[py, np.zeros(x.shape[1])]
        beta = np.linalg.lstsq(augmented, response, rcond=None)[0]
    gamma = np.linalg.lstsq(wz, wy-wx@beta, rcond=rank_tol)[0]
    fitted = x@beta + z@gamma
    residual = y-fitted
    wr = np.linalg.solve(lower, residual)
    return dict(beta=beta, nuisance=gamma, fitted=fitted, residual=residual,
                weighted_sse=float(wr@wr), xrank=xrank, zrank=zrank)


def reml_grid(matrix, y, rows, runs, model='arma11', gap='preserve'):
    """Independent QR-based restricted-likelihood grid reference (SciPy).

    Returns q x 6 [a,b,sigma2,score,rank,boundary] and grid scores. Columns of
    y are independent responses. Pilot SVD fixes a common orthonormal mean
    span; weighted QR supplies determinant and residuals. Zero variance and
    saturated pilots fail. Grid order matches the public deterministic contract.
    """
    from scipy.linalg import solve_triangular, cholesky, qr, svd
    matrix, y = np.asarray(matrix, float), np.asarray(y, float)
    u, s, _ = svd(matrix, full_matrices=False)
    rank = int(np.sum(s > 1e-10*s[0])) if len(s) else 0
    n = len(matrix)
    if rank >= n:
        raise ValueError('saturated pilot')
    basis = u[:, :rank]
    grid = [(a/10, b/10) for a in (range(-9, 10) if model != 'white' else [0])
            for b in (range(-9, 10) if model == 'arma11' else [0])]
    scores, variances = [], []
    for a, b in grid:
        r = correlation(rows, runs, model, a, b, gap)
        lower = cholesky(r, lower=True)
        weighted = solve_triangular(lower, basis, lower=True)
        q, upper = qr(weighted, mode='economic')
        wy = solve_triangular(lower, y, lower=True)
        residual = wy-q@(q.T@wy)
        variance = np.sum(residual**2, axis=0)/(n-rank)
        if np.any(variance <= np.finfo(float).tiny):
            raise ValueError('zero variance')
        determinant = 2*np.log(np.diag(lower)).sum()
        if rank:
            determinant += 2*np.log(np.abs(np.diag(upper))).sum()
        scores.append(determinant+(n-rank)*np.log(variance))
        variances.append(variance)
    scores, variances = np.array(scores), np.array(variances)
    index = np.argmin(scores, axis=0)
    result = np.zeros((y.shape[1], 6))
    for v, i in enumerate(index):
        a, b = grid[i]
        result[v] = [a, b, variances[i, v], scores[i, v], rank,
                     int(abs(a) == .9 or abs(b) == .9)]
    return result, scores
