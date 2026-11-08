#!/usr/bin/env python3
"""Crossed noise/fit benchmark with explicit train-only covariance estimation.

Default is the declared small pilot, not the precision-qualified full study.
--full uses the frozen dimensions; --adaptive expands cumulative replicates;
reports always state whether the MCSE target is met. No automatic claim of
scientific superiority. Requires rebuilt Picard/REMLfit, NumPy, SciPy, nibabel.
"""
import argparse
import hashlib
import json
import os
import resource
import sys
from pathlib import Path
import subprocess
import tempfile
import time
import nibabel as nib
import numpy as np
from scipy.linalg import cholesky, solve_triangular, svd
from scipy.optimize import brentq
from picard_covariance_reference import correlation, reml_grid


def write_matrix(path, matrix, rows, nfull, starts, roles):
    with open(path, 'w') as f:
        f.write(f'# <matrix\n# ni_type = "{matrix.shape[1]}*float"\n# ni_dimen = "{len(rows)}"\n')
        attrs = {'ColumnLabels': ';'.join(f'c{i}' for i in range(matrix.shape[1])),
                 'ColumnGroups': ','.join(str(1 if x else -1) for x in roles),
                 'GoodList': ','.join(map(str, rows)), 'NRowFull': nfull,
                 'RunStart': ','.join(map(str, starts)), 'RowTR': 1.5}
        for k, v in attrs.items(): f.write(f'# {k} = "{v}"\n')
        f.write('# >\n'); np.savetxt(f, matrix, fmt='%.9g'); f.write('# </matrix>\n')


def volume(path, data):
    nib.save(nib.Nifti1Image(data.T[:, None, None, :].astype(np.float32), np.eye(4)), path)


def bucket(path):
    img = nib.load(str(path)+'+orig.HEAD')
    return np.asarray(img.dataobj).reshape(img.shape[0], -1).T, img.header.info


def generator(config, smoke):
    """Fixed design/truth/censor seed streams; only noise varies by replicate."""
    seeds = np.random.SeedSequence(config['seed']).spawn(4)
    design_rng, beta_rng, censor_rng = [np.random.default_rng(s) for s in seeds[:3]]
    nr, nt, nc = (4, 64, 4) if smoke else (6, 240, 24)
    perrun = 2*nc; p = nr*perrun; n = nr*nt
    x, z = np.zeros((n, p)), np.zeros((n, 4*nr))
    condition = np.empty(p, int); trialrun = np.repeat(np.arange(nr), perrun)
    sample = np.arange(0, 32.0001, .1)
    hrf = sample**5*np.exp(-sample)/120
    hrf /= hrf.max()
    for r in range(nr):
        onsets = np.r_[12., 12+np.cumsum(design_rng.uniform(4, 8, perrun-1))]
        condition[r*perrun:(r+1)*perrun] = design_rng.permutation(np.tile(np.arange(nc), 2))
        for j, onset in enumerate(onsets):
            x[r*nt:(r+1)*nt, r*perrun+j] = np.interp(np.arange(nt)*1.5-onset, sample, hrf, left=0, right=0)
        z[r*nt:(r+1)*nt, r*4:(r+1)*4] = np.polynomial.legendre.legvander(np.linspace(-1, 1, nt), 3)
    means = beta_rng.normal(1, .5, nc)
    beta = means[condition]+beta_rng.normal(0, .5, p)
    compact = np.zeros((n, nr*nc))
    for j in range(p): compact[:, trialrun[j]*nc+condition[j]] += x[:, j]
    masks = {'none': np.arange(n)}
    remove = np.concatenate([r*nt+censor_rng.choice(nt, int(round(.1*nt)), replace=False) for r in range(nr)])
    masks['random'] = np.setdiff1d(np.arange(n), remove)
    blocked = []
    for r in range(nr):
        length = 4 if smoke else 12
        first = int(censor_rng.integers(0, nt-2*length))
        choices = [v for v in range(nt-length+1) if v+length <= first or v >= first+length]
        second = int(censor_rng.choice(choices))
        blocked.extend(r*nt+np.r_[np.arange(first, first+length), np.arange(second, second+length)])
    masks['blocks'] = np.setdiff1d(np.arange(n), blocked)
    # Match the AFNI float xmat boundary before producing synthetic signals.
    x, z, compact = [a.astype(np.float32).astype(float) for a in (x, z, compact)]
    return dict(x=x, z=z, compact=compact, beta=beta, condition=condition,
                trialrun=trialrun, nr=nr, nt=nt, masks=masks, noise_seed=seeds[3])


def draw_noise(g, truth, replicates):
    """Stationary independent runs; extending replicates preserves prior draws."""
    nr, nt = g['nr'], g['nt']
    rng = np.random.default_rng(g['noise_seed'])
    innovations = rng.normal(size=(replicates, nr*nt)).T
    covariance = correlation(np.arange(nt), np.zeros(nt, int), **truth)
    # Explicitly select the factor triangle: some LAPACK wrappers leave
    # scratch values in the unused triangle even for cholesky(), and dense
    # multiplication (unlike triangular solves) would consume those values.
    lower = np.tril(cholesky(covariance, lower=True))
    return np.vstack([lower@innovations[r*nt:(r+1)*nt] for r in range(nr)])


def fit_path(x, z, y, rows, runs, parameters, model, gap, fractions):
    """Independent SVD fractional GLS path, grouped by exact covariance pair."""
    p, nv = x.shape[1], y.shape[1]
    result = np.zeros((len(fractions), p, nv))
    active = np.flatnonzero(np.any(x != 0, axis=0))
    nuisance = z[:, np.any(z != 0, axis=0)]
    for pair in np.unique(parameters, axis=0):
        vox = np.flatnonzero(np.all(parameters == pair, axis=1))
        lower = cholesky(correlation(rows, runs, model, *pair, gap), lower=True)
        wx = solve_triangular(lower, x[:, active], lower=True)
        wy = solve_triangular(lower, y[:, vox], lower=True)
        wz = solve_triangular(lower, nuisance, lower=True)
        q = svd(wz, full_matrices=False)[0] if wz.shape[1] else np.empty((len(rows), 0))
        px, py = wx-q@(q.T@wx), wy-q@(q.T@wy)
        u, s, vt = svd(px, full_matrices=False)
        keep = s > 1e-10*s[0]
        u, s, vt = u[:, keep], s[keep], vt[keep]
        rotation = u.T@py
        reference = rotation/s[:, None]
        for k, fraction in enumerate(fractions):
            for c, v in enumerate(vox):
                if fraction == 1: coefficients = reference[:, c]
                else:
                    target = fraction*np.linalg.norm(reference[:, c])
                    if target == 0: coefficients = np.zeros(len(s))
                    else:
                        def difference(alpha):
                            return np.linalg.norm(s*rotation[:, c]/(s*s+alpha))-target
                        high = 1.
                        while difference(high) > 0: high *= 2
                        alpha = brentq(difference, 0, high, xtol=1e-12)
                        coefficients = s*rotation[:, c]/(s*s+alpha)
                result[k, active, v] = vt.T@coefficients
    return result


def paired_summary(values):
    values = np.asarray(values, float)
    return {'mean': float(values.mean()),
            'mcse': float(values.std(ddof=1)/np.sqrt(len(values))) if len(values)>1 else None,
            'n': len(values)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--picard', required=True)
    parser.add_argument('--remlfit', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--estimator', choices=['reference', 'native'], default='reference')
    parser.add_argument('--full', action='store_true')
    parser.add_argument('--adaptive', action='store_true',
                        help='full protocol: expand cumulative replicates to MCSE target or cap')
    parser.add_argument('--replicates', type=int)
    parser.add_argument('--censor', choices=['none', 'random', 'blocks', 'all'], default='random')
    args = parser.parse_args()
    if args.replicates is None: args.replicates = 200 if args.full else 8
    if args.replicates < 2: parser.error('at least two replicates required')
    picard, reml = str(Path(args.picard).resolve()), str(Path(args.remlfit).resolve())
    protocol = Path(__file__).with_name('picard_ng1_benchmark.json')
    protocol_bytes = protocol.read_bytes()
    config = json.loads(protocol_bytes)
    if args.adaptive:
        if not args.full or args.censor != 'all':
            parser.error('--adaptive requires --full --censor all')
        if args.replicates < config['replicates_initial']:
            parser.error('adaptive full study starts with at least 200 replicates')
        if args.replicates > config['replicates_max']:
            parser.error('replicate count exceeds the frozen adaptive cap')
        history = []
        adaptive_start = time.perf_counter()
        counts = list(range(args.replicates, config['replicates_max'],
                            config['replicate_batch']))+[config['replicates_max']]
        for count in counts:
            output = Path(args.output).with_suffix(f'.reps_{count}.json')
            # A failed child must never be mistaken for an earlier report.
            output.unlink(missing_ok=True)
            command = [sys.executable, str(Path(__file__).resolve()),
                       '--picard', picard, '--remlfit', reml,
                       '--estimator', args.estimator, '--full', '--censor', 'all',
                       '--replicates', str(count), '--output', str(output)]
            completed = subprocess.run(command, check=False)
            if not output.exists():
                raise RuntimeError(f'adaptive batch failed without report: {completed.returncode}')
            report = json.loads(output.read_text())
            if completed.returncode and not report['failures']:
                raise RuntimeError(f'adaptive child exited {completed.returncode}')
            history.append({'replicates': count, 'report': str(output),
                            'precision_target_met': report['precision_target_met'],
                            'failed_cells': len(report['failures'])})
            if completed.returncode or report['precision_target_met']:
                break
        report['adaptive_history'] = history
        report['adaptive_total_wall_seconds'] = time.perf_counter()-adaptive_start
        report['adaptive_note'] = ('Cumulative reruns preserve earlier realizations; '
                                   'precision stopping is not an inferential significance test.')
        Path(args.output).write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        if report['failures']: raise SystemExit(1)
        return
    g = generator(config, not args.full)
    nr, nt = g['nr'], g['nt']; n = nr*nt
    fractions = config['ridge_fraction_grid'] if args.full else [1., .7, .4]
    truths = config['noise_cases'] if args.full else [config['noise_cases'][i] for i in [0, 2, 4]]
    tuning_runs = np.arange(nr-2 if args.full else nr-1)
    evaluation_runs = np.arange(len(tuning_runs), nr)
    env = dict(os.environ, AFNI_FLOATSCAN='NO', AFNI_WRITE_1D_AS_PREFIX='YES',
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    records, failures, checks, costs = [], [], [], []
    begin = time.perf_counter()
    command_counter = 0
    with tempfile.TemporaryDirectory(prefix='picard_benchmark_') as tmp:
        folder = Path(tmp)
        def execute(binary, options):
            nonlocal command_counter
            command_counter += 1
            started = time.perf_counter()
            command = [binary, *map(str, options)]
            p = subprocess.Popen(command, env=env, cwd=folder, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            # communicate drains output; per-call time is recorded. Peak RSS is
            # measured for the whole harness externally, not guessed per cell.
            out, err = p.communicate()
            costs.append({'program': Path(binary).name, 'seconds': time.perf_counter()-started,
                          'returncode': p.returncode})
            if p.returncode: raise RuntimeError(err[-2500:])
            return out
        for truth_id, truth in enumerate(truths):
            # Same innovations across truth families; exact stationary full-run covariance.
            noise = draw_noise(g, truth, args.replicates)
            yfull = (g['x']@g['beta'])[:, None]+noise
            yfull = yfull.astype(np.float32).astype(float)
            volume(folder/'data.nii', yfull)
            for censor, rows in g['masks'].items():
                if args.censor != 'all' and censor != args.censor: continue
                runs = rows//nt; x, z, y = g['x'][rows], g['z'][rows], yfull[rows]
                train = np.isin(runs, tuning_runs); test = np.isin(runs, evaluation_runs)
                test_trials = np.isin(g['trialrun'], evaluation_runs)
                for gap in ['preserve', 'reset']:
                    for pilot_name in ['full', 'compact']:
                        pilot = np.c_[x if pilot_name == 'full' else g['compact'][rows], z]
                        pilot_file = folder/'pilot.1D'
                        write_matrix(pilot_file, pilot, rows, n, np.arange(nr)*nt, np.ones(pilot.shape[1], int))
                        simple = np.c_[np.ones(len(rows)), np.linspace(-1, 1, len(rows))]
                        write_matrix(folder/'simple.1D', simple, rows, n, np.arange(nr)*nt, [0, 1])
                        for model in ['white', 'ar1', 'arma11']:
                            cell = dict(truth=truth, truth_id=truth_id, fitted=model, gap=gap,
                                        censor=censor, pilot=pilot_name, estimator=args.estimator)
                            try:
                                def estimate(selected, external=False):
                                    if external:
                                        path = folder/f'external{command_counter}'
                                        active = np.any(pilot[selected] != 0, axis=0)
                                        write_matrix(folder/'external.1D', pilot[selected][:, active], rows[selected], n,
                                                     np.arange(nr)*nt, np.ones(sum(active), int))
                                        execute(reml, ['-input', folder/'data.nii', '-matrix', folder/'external.1D',
                                                       '-NEGcor', '-Grid', 3, '-MAXa', .9, '-MAXb', .9,
                                                       '-CORcut', 0, '-Rvar', path])
                                        return bucket(path)[0][:2].T
                                    if args.estimator == 'reference':
                                        return reml_grid(pilot[selected], y[selected], rows[selected], runs[selected], model, gap)[0][:, :2]
                                    np.savetxt(folder/'rows.1D', rows[selected], fmt='%d')
                                    path = folder/f'native{command_counter}'
                                    execute(picard, ['-input', folder/'data.nii', '-matrix', folder/'simple.1D',
                                        '-covariance', model, '-covariance_gaps', gap,
                                        '-covariance_estimate', pilot_file, '-covariance_rows', folder/'rows.1D',
                                        '-fraction', 1, '-prefix', path])
                                    return bucket(path)[0][-2:].T
                                losses = np.zeros((len(fractions), args.replicates)); pairs = 0
                                for held in tuning_runs:
                                    tr = train & (runs != held); va = runs == held
                                    pars = estimate(tr)
                                    trained = fit_path(x[tr], z[tr], y[tr], rows[tr], runs[tr], pars, model, gap, fractions)
                                    reference = fit_path(x[va], z[va], y[va], rows[va], runs[va], pars, model, gap, [1])[0]
                                    for j in np.flatnonzero(np.isin(g['trialrun'], tuning_runs) & (g['trialrun'] != held)):
                                        matches = np.flatnonzero((g['trialrun'] == held) & (g['condition'] == g['condition'][j]))
                                        losses += ((trained[:, j, None, :]-reference[matches][None])**2).sum(axis=1)
                                        pairs += len(matches)
                                losses /= pairs
                                chosen = np.argmin(losses, axis=0)  # first/largest fraction wins ties
                                pars = estimate(train)
                                fitted = fit_path(x[test], z[test], y[test], rows[test], runs[test], pars, model, gap, fractions)
                                picked = fitted[chosen, :, np.arange(args.replicates)].T
                                # Exercise actual Picard fitting on evaluation rows using frozen parameters.
                                active_x = np.flatnonzero(np.any(x[test] != 0, axis=0))
                                active_z = np.flatnonzero(np.any(z[test] != 0, axis=0))
                                evalmat = np.c_[x[test][:, active_x], z[test][:, active_z]]
                                write_matrix(folder/'eval.1D', evalmat, rows[test], n, np.arange(nr)*nt,
                                             [1]*len(active_x)+[0]*len(active_z))
                                volume(folder/'parameters.nii', pars.T)
                                for k in np.unique(chosen):
                                    path = folder/f'eval{command_counter}'
                                    options = ['-input', folder/'data.nii', '-matrix', folder/'eval.1D',
                                               '-covariance', model, '-covariance_gaps', gap,
                                               '-fraction', fractions[k], '-prefix', path]
                                    if model != 'white': options += ['-covariance_map', folder/'parameters.nii']
                                    execute(picard, options)
                                    actual = bucket(path)[0][:len(active_x)]
                                    np.testing.assert_allclose(actual, fitted[k, active_x], atol=2e-4, rtol=2e-4)
                                checks.append(cell)
                                for name, betas in [('unpenalized', fitted[0]), ('tuned_ridge', picked)]:
                                    estimate_beta = betas[test_trials]
                                    error = estimate_beta-g['beta'][test_trials, None]
                                    mse = np.mean(error**2, axis=0)
                                    mean_error = error.mean(axis=1)
                                    target = g['beta'][test_trials]
                                    centered = target-target.mean()
                                    slope = centered@estimate_beta/(centered@centered)
                                    intercept = estimate_beta.mean(axis=0)-slope*target.mean()
                                    true_r = correlation(rows[test], runs[test], **truth)
                                    covariance_error = [float(np.linalg.norm(
                                        correlation(rows[test], runs[test], model, *pair, gap)-true_r)
                                        /np.linalg.norm(true_r)) for pair in pars]
                                    records.append(dict(cell, method=name, mse=mse.tolist(),
                                        amplitude_slope=slope.tolist(), amplitude_intercept=intercept.tolist(),
                                        relative_covariance_error=covariance_error,
                                        summary=paired_summary(mse), conditional_bias2=float(np.mean(mean_error**2)),
                                        conditional_variance=float(np.mean(np.var(estimate_beta, axis=1, ddof=0))),
                                        selected_fraction=[fractions[i] for i in chosen], parameters=pars.tolist(),
                                        boundary_count=int(np.sum(np.any(np.abs(pars) >= .89999, axis=1)))))
                                # Matched true-covariance oracle is separate from deployable methods.
                                if model == truth['model'] and gap == 'preserve':
                                    oracle = np.tile([truth.get('ar', 0), truth.get('ma', 0)], (args.replicates, 1))
                                    oracle_fit = fit_path(x[test], z[test], y[test], rows[test], runs[test], oracle, model, gap, [1])[0]
                                    records.append(dict(cell, method='known_covariance_GLS',
                                        mse=np.mean((oracle_fit[test_trials]-g['beta'][test_trials, None])**2, axis=0).tolist()))
                                # AFNI reference uses its own grid, preserve-gap policy only.
                                if model == 'arma11' and gap == 'preserve':
                                    ep = estimate(train, external=True)
                                    external_fit = fit_path(x[test], z[test], y[test], rows[test], runs[test], ep, model, gap, [1])[0]
                                    records.append(dict(cell, method='AFNI_REML_GLS', parameters=ep.tolist(),
                                        mse=np.mean((external_fit[test_trials]-g['beta'][test_trials, None])**2, axis=0).tolist()))
                            except Exception as exc:
                                failures.append(dict(cell, error=str(exc)))
        contrasts = []
        for record in records:
            if record['method'] != 'tuned_ridge': continue
            keys = ['truth_id', 'fitted', 'gap', 'censor', 'pilot', 'estimator']
            baseline = next(r for r in records if r['method'] == 'unpenalized' and all(r[k] == record[k] for k in keys))
            contrasts.append(dict({k: record[k] for k in keys}, contrast='ridge_minus_GLS',
                                  **paired_summary(np.array(record['mse'])-baseline['mse'])))
        # Paired comparisons retain exactly the same design, truth, noise and
        # censor realization; never pair results from different pilots.
        for record in records:
            if record['method'] not in ('unpenalized', 'tuned_ridge'): continue
            for kind, field, current, baseline_value in [
                    ('ARMA_minus_AR1', 'fitted', 'arma11', 'ar1'),
                    ('preserve_minus_reset', 'gap', 'preserve', 'reset')]:
                if record[field] != current: continue
                keys = ['truth_id', 'fitted', 'gap', 'censor', 'pilot', 'estimator', 'method']
                match = next((r for r in records if r[field] == baseline_value
                              and all(r[k] == record[k] for k in keys if k != field)), None)
                if match:
                    contrasts.append(dict({k: record[k] for k in keys}, contrast=kind,
                        **paired_summary(np.array(record['mse'])-match['mse'])))
    result = dict(schema_version=1, protocol_sha256=hashlib.sha256(protocol_bytes).hexdigest(),
                  profile='full_dimensions' if args.full else 'reduced_pilot', estimator=args.estimator,
                  deviations=[] if args.full else ['4 runs x 64 TRs', '4 conditions twice/run',
                    '3 tuning runs and 1 evaluation run', '3 fractions', '3 noise cases',
                    'censor patterns selected by --censor; fixed replicate count, no adaptive expansion'],
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  git_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=Path(__file__).resolve().parents[3], text=True).strip(),
                  max_child_rss_raw=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                  rss_units='bytes' if sys.platform == 'darwin' else 'KiB',
                  thread_settings={'OMP_NUM_THREADS': 1, 'OPENBLAS_NUM_THREADS': 1},
                  replicates=args.replicates, censor=args.censor, tuning_runs=tuning_runs.tolist(),
                  evaluation_runs=evaluation_runs.tolist(), records=records, contrasts=contrasts,
                  failures=failures, validated_picard_cells=len(checks), command_costs=costs,
                  wall_seconds=time.perf_counter()-begin, executables=[picard, reml],
                  executable_sha256=[hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in [picard, reml]],
                  precision_target_met=bool(contrasts) and not failures and all(c['mcse'] <= .01 for c in contrasts),
                  precision_scope='paired contrasts in executed cells only; reduced runs do not certify full protocol',
                  full_protocol_dimensions=args.full,
                  full_protocol_coverage=args.full and args.censor == 'all' and args.replicates >= 200,
                  afni_reference='NEGcor Grid3 MAXa=.9 MAXb=.9 CORcut=0; ARMA only',
                  note='Pilot evidence, not superiority or completed full scientific study.')
    target = Path(args.output)
    target.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(f'{len(records)} results, {len(checks)} Picard-validated cells, {len(failures)} failures; {target}')
    if failures: raise SystemExit(1)


if __name__ == '__main__': main()
