#!/usr/bin/env python3
"""Small NG2 ordering/ablation probe, not an independent prediction benchmark.

Usage: python benchmark_picard_ng2.py PICARD DECONVOLVE --output REPORT.json
The second binary argument follows the existing fixture-helper convention;
it is not used to estimate coefficients. Only Picard is executed in this probe.

Design, trial truth and HRF identities stay fixed across noise realizations.
Compare pilot-covariance -> HRF -> PCs -> ridge with a limited alternative:
white-fit HRF -> selected-HRF pilot covariance -> PCs -> ridge. Both reuse
full data for learning, so neither supplies an independent generalization score.
The alternative is a reference-only research probe, not a new production mode.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path

# This companion reuses the independent NumPy/augmented-system reference and
# fixture writers from the integration test. Importing does not run its tests.
from test_picard_ng2 import (np, os, tempfile, BINARY, candidate_manifest,
    xmat, volume, run, bucket, correlation, reml_grid, pca_reference, staged_reference)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('picard')
    parser.add_argument('deconvolve')
    parser.add_argument('--output', required=True)
    parser.add_argument('--replicates', type=int, default=8)
    parser.add_argument('--model', choices=['ar1', 'arma11'], default='ar1')
    args = parser.parse_args()
    if args.replicates < 2: parser.error('at least two paired realizations required')
    destination = Path(args.output).resolve()
    start = time.perf_counter()
    rng = np.random.default_rng(20261109)
    nr, nt, ns, nv = 3, 36, 4, 10
    n = nr*nt; p = 3*nr
    rows = np.arange(n); runs = rows//nt
    trial_runs = np.repeat(np.arange(nr), 3)
    conditions = np.tile(np.arange(3), nr)
    roles = [1]*p+[0]*nr
    design = []
    for delay, width in [(0., 2.), (1.5, 3.5)]:
        matrix = np.zeros((n, p+nr))
        for r in range(nr):
            for trial, onset in enumerate([5, 15, 25]):
                matrix[r*nt:(r+1)*nt, r*3+trial] = np.exp(-.5*((np.arange(nt)-onset-delay)/width)**2)
            matrix[r*nt:(r+1)*nt, p+r] = 1
        design.append(matrix)
    true_hrf = np.arange(ns) % 2
    truth = np.tile([1., 2., -1.], nr)[:, None]+rng.normal(0, .25, (p, ns))
    loadings = rng.normal(size=(2, nv))*.5
    true_covariance = correlation(rows, runs, 'arma11', .6, .2)
    lower = np.linalg.cholesky(true_covariance)
    records, failures = [], []
    with tempfile.TemporaryDirectory(prefix='picard-ng2-order-') as work:
        os.chdir(work)
        manifest, matrices = candidate_manifest('hrf', design, roles, rows,
            np.arange(nr)*nt, ids=(11,4), trial_runs=trial_runs, conditions=conditions)
        pilot = xmat('pilot.1D', matrices[0], roles, rows, n, np.arange(nr)*nt)
        volume('pool.nii', [np.arange(nv) >= ns])
        volume('analysis.nii', [np.arange(nv) < ns])
        np.savetxt('trials.1D', np.c_[np.arange(p),np.arange(p),conditions,trial_runs], fmt='%d')
        for replicate in range(args.replicates):
            # Each realization gets its own noise seed; adding replicates
            # leaves the complete earlier realization (including PCs) intact.
            noise_rng = np.random.default_rng(np.random.SeedSequence([20261110, replicate]))
            response = lower@noise_rng.normal(size=(n, nv))*.5
            response += (lower@noise_rng.normal(size=(n, 2)))@loadings
            for voxel, hrf in enumerate(true_hrf):
                response[:, voxel] += matrices[hrf][:, :p]@truth[:, voxel]
            response = response.astype(np.float32).astype(float)
            try:
                pcs = pca_reference(response, runs, np.arange(ns, nv), 2)
                parameters = reml_grid(pilot, response[:, :ns], rows, runs, args.model)[0][:, :2]
                common = dict(matrices=matrices, response=response[:, :ns], roles=roles,
                    good=rows, runs=runs, gap='preserve', pcs=pcs, trial_runs=trial_runs,
                    conditions=conditions, select_voxels=np.arange(ns))
                primary = staged_reference(**common, parameters=parameters, model=args.model, fractions=[1,.7,.3])
                white = staged_reference(**common, parameters=np.zeros((ns,2)), model='white', fractions=[1,.7,.3])
                alternative_parameters = np.zeros((ns, 2))
                # Estimate only AFTER locking each voxel's white-fit HRF.
                # REML scores are never compared across those mean-model spans.
                for hrf in range(2):
                    voxels = np.flatnonzero(white['hrf'] == hrf)
                    if len(voxels):
                        alternative_parameters[voxels] = reml_grid(matrices[hrf],
                            response[:, voxels], rows, runs, args.model)[0][:, :2]
                alternative = staged_reference(**common, parameters=alternative_parameters,
                    model=args.model, fractions=[1,.7,.3], locked_hrf=white['hrf'])
                no_pc = staged_reference(**common, parameters=parameters, model=args.model,
                                         fractions=[1,.7,.3], fixed_count=0)
                no_ridge = staged_reference(**common, parameters=parameters, model=args.model,
                                            fractions=[1], fixed_count=primary['chosen'])
                # Check the actual new CLI for every realization of the primary
                # ordering. The other paths are explicitly labeled references.
                volume('data.nii', response)
                prefix = f'fit{replicate}'
                run([BINARY,'-input','data.nii','-mask','analysis.nii',
                    '-candidate_manifest',manifest,'-candidate_policy','hrf',
                    '-single_trial','-trial_info','trials.1D','-noise_pool','pool.nii',
                    '-pc_counts','0,1,2','-pc_mask','analysis.nii','-fractions','1,.7,.3',
                    '-covariance',args.model,'-covariance_estimate','pilot.1D',
                    '-covariance_freeze','-prefix',prefix])
                actual, header, labels = bucket(prefix)
                np.testing.assert_allclose(actual[:p,:ns], primary['beta'], rtol=2e-4, atol=2e-4)
                assert header['PICARD_PC_COUNT'] == primary['chosen']
                np.testing.assert_array_equal(actual[labels['Picard_DesignIndex'],:ns], primary['hrf'])
                for name, result in [('covariance_first',primary),('HRF_first_reference',alternative),
                                     ('white_reference',white),('no_PC_reference',no_pc),
                                     ('unpenalized_reference',no_ridge)]:
                    records.append(dict(replicate=replicate, method=name,
                        beta_mse=float(np.mean((result['beta']-truth)**2)),
                        HRF_error=float(np.mean(result['hrf'] != true_hrf)),
                        PC_count=result['chosen'], selected_fraction_index=list(map(int,result['wins']))))
            except Exception as exc:
                failures.append(dict(replicate=replicate, error=str(exc)))
    summaries = {}
    for method in sorted({record['method'] for record in records}):
        values = [record['beta_mse'] for record in records if record['method'] == method]
        summaries[method] = dict(mean_beta_mse=float(np.mean(values)), n=len(values))
    paired = []
    for replicate in range(args.replicates):
        match = {record['method']: record['beta_mse'] for record in records if record['replicate'] == replicate}
        if len(match) == 5:
            paired.append(match['covariance_first']-match['HRF_first_reference'])
    report = dict(schema_version=1, seed_design=20261109, seed_noise=20261110,
        replicates=args.replicates, fitted_covariance=args.model,
        truth_covariance=dict(model='arma11', ar=.6, ma=.2),
        profile='3 runs x 36 TR, 3 conditions/run, 4 signal and 6 pool voxels; no censoring',
        selection='full-data conditional adaptive fits; no independent test partition',
        purpose='limited ordering/ablation probe; not an optimal-order or superiority claim',
        records=records, failures=failures, summaries=summaries,
        paired_covariance_first_minus_HRF_first=dict(n=len(paired),
            mean=float(np.mean(paired)) if paired else None,
            mcse=float(np.std(paired,ddof=1)/np.sqrt(len(paired))) if len(paired)>1 else None),
        wall_seconds=time.perf_counter()-start,
        executable_sha256=hashlib.sha256(Path(BINARY).read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    destination.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(f'{len(paired)} validated realizations; {len(failures)} failures; {destination}')
    if failures: raise SystemExit(1)


if __name__ == '__main__':
    main()
