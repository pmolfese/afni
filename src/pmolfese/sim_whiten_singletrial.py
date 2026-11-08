"""
Toy pilot: does noise-model-aware (whitened) ridge improve single-trial beta
recovery over white-noise ridge, and does it matter HOW the AR parameter is
estimated when the single-trial design is nearly saturated?

Design: 4 runs x 240 TRs (TR=1.5 s), jittered ISI, one regressor per trial
(~80/run), per-run Legendre drift (order 3, unpenalized), AR(1) or ARMA(1,1)
noise, true trial betas ~ N(1, 0.5^2), noise sd = 1.

Estimators (all ridge on trial columns, drift unpenalized, lambda chosen as the
single best value averaged over realizations, so tuning is not a confound):
  OLS-white     : no whitening
  GLS-oracle    : whiten with the TRUE noise covariance
  GLS-estSTD    : AR(1) estimated by REML (profile over a) from residuals of
                  the full single-trial model (what a naive REML-on-X would do)
  GLS-estONOFF  : AR(1) estimated by REML from a compact condition-level
                  (ONOFF + drift) model
Not a validation of anything in 3dPicard -- an effect-size / bias probe.
"""
import numpy as np
import math


def lfilter(b, a, x):
    """numpy-only IIR filter, a[0]==1 (scipy is unavailable in this sandbox)."""
    y = np.zeros(len(x))
    for n in range(len(x)):
        acc = sum(b[k] * x[n - k] for k in range(len(b)) if n - k >= 0)
        acc -= sum(a[k] * y[n - k] for k in range(1, len(a)) if n - k >= 0)
        y[n] = acc
    return y


def toeplitz(c):
    n = len(c)
    i = np.arange(n)
    return np.asarray(c)[np.abs(i[:, None] - i[None, :])]


def cholesky(A, lower=True):
    return np.linalg.cholesky(A)


def solve_triangular(L, B, lower=True):
    return np.linalg.solve(L, B)


class gam:
    @staticmethod
    def pdf(t, k):
        t = np.asarray(t, float)
        return np.where(t > 0, t ** (k - 1) * np.exp(-t) / math.gamma(k), 0.0)

rng = np.random.default_rng(7)
TR, NRUN, NT = 1.5, 4, 240
LAMS = np.logspace(-1, 3, 17)


def hrf(tr, dur=32.0):
    t = np.arange(0, dur, tr)
    h = gam.pdf(t, 6) - gam.pdf(t, 16) / 6.0
    return h / np.abs(h).max()


def make_design(rng):
    cols, runs_of_col = [], []
    for r in range(NRUN):
        t = 2.0
        onsets = []
        while t < NT * TR - 20:
            onsets.append(t)
            t += rng.uniform(3.0, 6.0)
        h = hrf(TR)
        for o in onsets:
            sig = np.zeros(NRUN * NT)
            # onset on TR grid with sub-TR linear split
            idx = o / TR
            i0 = int(np.floor(idx)); fr = idx - i0
            for k, w in ((i0, 1 - fr), (i0 + 1, fr)):
                seg = h[: max(0, NT - k)] * w
                sig[r * NT + k: r * NT + k + len(seg)] += seg
            cols.append(sig); runs_of_col.append(r)
    X = np.array(cols).T
    Z = np.zeros((NRUN * NT, NRUN * 4))
    tt = np.linspace(-1, 1, NT)
    for r in range(NRUN):
        for d in range(4):
            Z[r * NT:(r + 1) * NT, r * 4 + d] = np.polynomial.legendre.Legendre.basis(d)(tt)
    return X, Z, np.array(runs_of_col)


def ar1_whiten(M, a):
    """Per-run AR(1) whitening; M is (NRUN*NT, k)."""
    out = np.empty_like(M)
    for r in range(NRUN):
        s = slice(r * NT, (r + 1) * NT)
        v = M[s]
        w = np.empty_like(v)
        w[0] = np.sqrt(1 - a * a) * v[0]
        w[1:] = v[1:] - a * v[:-1]
        out[s] = w
    return out


def arma_cov_chol(a, b):
    """Cholesky factor of per-run ARMA(1,1) covariance (unit innovations)."""
    imp = lfilter([1, b], [1, -a], np.r_[1.0, np.zeros(400)])
    acv = np.array([imp[:len(imp) - k] @ imp[k:] for k in range(NT)])
    L = cholesky(toeplitz(acv), lower=True)
    return L


def gen_noise(kind, par, rng):
    e = np.empty(NRUN * NT)
    if kind == "ar1":
        a = par
        for r in range(NRUN):
            x = lfilter([1], [1, -a], rng.standard_normal(NT + 200))[200:]
            e[r * NT:(r + 1) * NT] = x
    else:
        a, b = par
        for r in range(NRUN):
            x = lfilter([1, b], [1, -a], rng.standard_normal(NT + 200))[200:]
            e[r * NT:(r + 1) * NT] = x
    return e / e.std()


def whiten_oracle_fn(kind, par):
    if kind == "ar1":
        return lambda M: ar1_whiten(M, par)
    L = arma_cov_chol(*par)
    def f(M):
        out = np.empty_like(M)
        for r in range(NRUN):
            s = slice(r * NT, (r + 1) * NT)
            out[s] = solve_triangular(L, M[s], lower=True)
        return out
    return f


def reml_ar1(y, D):
    """Profile REML over AR(1) coefficient a for response y, design D."""
    best, besta = -np.inf, 0.0
    n, p = D.shape
    for a in np.linspace(-0.3, 0.9, 13):
        Dw, yw = ar1_whiten(D, a), ar1_whiten(y[:, None], a)[:, 0]
        coef, *_ = np.linalg.lstsq(Dw, yw, rcond=None)
        rss = np.sum((yw - Dw @ coef) ** 2)
        sign, ld = np.linalg.slogdet(Dw.T @ Dw)
        ll = NRUN * 0.5 * np.log(1 - a * a) - 0.5 * (n - p) * np.log(rss / (n - p)) - 0.5 * ld
        if ll > best:
            best, besta = ll, a
    return besta


def ridge_betas(Xw, Zw, yw, lams):
    """Project Zw out of Xw,yw, then ridge via SVD; returns (len(lams), p)."""
    Q, _ = np.linalg.qr(Zw)
    Xp = Xw - Q @ (Q.T @ Xw)
    yp = yw - Q @ (Q.T @ yw)
    U, s, Vt = np.linalg.svd(Xp, full_matrices=False)
    uty = U.T @ yp
    return np.array([Vt.T @ (s / (s * s + l) * uty) for l in lams])


def run_scenario(kind, par, nrep=80):
    names = ["OLS-white", "GLS-oracle", "GLS-estSTD", "GLS-estONOFF"]
    rmse = {k: np.zeros(len(LAMS)) for k in names}
    corr = {k: np.zeros(len(LAMS)) for k in names}
    ahat = {"GLS-estSTD": [], "GLS-estONOFF": []}
    X, Z, runs = make_design(rng)
    onoff = np.zeros((NRUN * NT, 1))
    onoff[:, 0] = X.sum(axis=1)
    orc = whiten_oracle_fn(kind, par)
    Xo, Zo = orc(X), orc(Z)
    for _ in range(nrep):
        beta = rng.normal(1.0, 0.5, X.shape[1])
        e = gen_noise(kind, par, rng)
        y = X @ beta + e
        # pre-estimate AR(1) from the two residual models (drift always in)
        a_std = reml_ar1(y, np.hstack([X, Z]))
        a_onoff = reml_ar1(y, np.hstack([onoff, Z]))
        ahat["GLS-estSTD"].append(a_std); ahat["GLS-estONOFF"].append(a_onoff)
        variants = {
            "OLS-white": (X, Z, y),
            "GLS-oracle": (Xo, Zo, orc(y[:, None])[:, 0]),
            "GLS-estSTD": (ar1_whiten(X, a_std), ar1_whiten(Z, a_std), ar1_whiten(y[:, None], a_std)[:, 0]),
            "GLS-estONOFF": (ar1_whiten(X, a_onoff), ar1_whiten(Z, a_onoff), ar1_whiten(y[:, None], a_onoff)[:, 0]),
        }
        for k, (Xw, Zw, yw) in variants.items():
            B = ridge_betas(Xw, Zw, yw, LAMS)
            rmse[k] += np.sqrt(((B - beta) ** 2).mean(axis=1))
            bc = B - B.mean(axis=1, keepdims=True)
            tc = beta - beta.mean()
            corr[k] += (bc @ tc) / (np.linalg.norm(bc, axis=1) * np.linalg.norm(tc))
    print(f"\n=== noise {kind} {par}: p={X.shape[1]} trials, n={NRUN*NT} rows ===")
    for k in names:
        r = rmse[k] / nrep; c = corr[k] / nrep
        extra = ""
        if k in ahat:
            extra = f"  mean a_hat={np.mean(ahat[k]):.3f} (sd {np.std(ahat[k]):.3f})"
        print(f"{k:13s} best-lambda RMSE={r.min():.4f} (lam={LAMS[r.argmin()]:.2g})"
              f"  best corr(beta_hat,beta)={c.max():.4f}{extra}")
    base = (rmse['OLS-white'] / nrep).min()
    for k in names[1:]:
        print(f"   RMSE change vs OLS-white, {k}: {100*((rmse[k]/nrep).min()/base-1):+.1f}%")


if __name__ == "__main__":
    import sys; sc = sys.argv[1]
    {"a": lambda: run_scenario("ar1", 0.3), "b": lambda: run_scenario("ar1", 0.6), "c": lambda: run_scenario("arma", (0.8, -0.5))}[sc]()
