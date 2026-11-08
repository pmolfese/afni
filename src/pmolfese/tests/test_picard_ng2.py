#!/usr/bin/env python3
"""NG2: independent staged GLS, shared PCA/count, and frozen-pilot integration.

The reference uses NumPy least squares/Cholesky and an augmented-system Brent
fraction solver. Only fixture I/O and existing independent oracles are reused;
no production Picard numerical function computes an expected answer here.
"""
from pathlib import Path
source = Path(__file__).with_name('test_picard.py').read_text()
exec(compile(source.split("with tempfile.TemporaryDirectory(prefix='picard-tests-')")[0],
             str(Path(__file__).with_name('test_picard.py')), 'exec'))
from picard_covariance_reference import correlation, reml_grid


def bucket(prefix):
    """Read a bucket by labels; native diagnostics change its total size."""
    img = nib.load(str(prefix)+'+orig.HEAD')
    data = np.asarray(img.dataobj).reshape(img.shape[0], -1).T
    info = img.header.info
    labels = {name: index for index, name in enumerate(info['BRICK_LABS'].split('~'))}
    return data, info, labels


def gls(matrix, response, roles, lower, fraction=1):
    """Independent fit in weighted space, reconstructed in original units."""
    wx = np.linalg.solve(lower, matrix)
    wy = np.linalg.solve(lower, response)
    beta, gamma, _, residual, _, _ = oracle(wx, wy, roles, fraction)
    mask = np.asarray(roles, bool)
    fitted = matrix[:, mask]@beta + matrix[:, ~mask]@gamma
    return beta, gamma, fitted, response-fitted, np.mean(residual**2, axis=0)


def pca_reference(data, runs, pool_voxels, count):
    """Original-time PCA by SVD of normalized responses, not a Gram SVD."""
    pcs = np.zeros((len(data), count))
    for run_id in np.unique(runs):
        local = data[runs == run_id][:, pool_voxels].copy()
        local -= local.mean(axis=0)
        local /= np.linalg.norm(local, axis=0)
        u = np.linalg.svd(local, full_matrices=False)[0][:, :count]
        u /= u.std(axis=0, ddof=1)
        u *= np.sign(u[np.argmax(np.abs(u), axis=0), np.arange(count)])
        pcs[runs == run_id] = u
    return pcs


def augmented_design(matrix, pcs, runs, count):
    """Append run-local nuisance PCs without changing task-column identities."""
    extra = np.zeros((len(matrix), len(np.unique(runs))*count))
    for ordinal, run_id in enumerate(np.unique(runs)):
        extra[runs == run_id, ordinal*count:(ordinal+1)*count] = pcs[runs == run_id, :count]
    return np.c_[matrix, extra]


def staged_reference(matrices, response, roles, good, runs, parameters,
                     model, gap, pcs, trial_runs, conditions, fractions,
                     select_voxels, fixed_count=None, locked_hrf=None):
    """Full conditional pipeline; optional locks support ordering probes.

    Noise is fixed per voxel across HRFs. PC losses are reduced across the
    complete selection mask only after every voxel has its own weighted fits.
    In leave-one-run-out tuning, block-diagonal whitening commutes with removing
    an entire run; this reference therefore uses the existing repeat oracle on
    independently whitened inputs. That shortcut is NOT valid for arbitrary
    within-run prediction folds.
    """
    nv = response.shape[1]
    ntrial = sum(roles)
    maxpc = pcs.shape[1]
    factors = [np.linalg.cholesky(correlation(good, runs, model, *pair, gap))
               for pair in parameters]
    scores = np.empty((len(matrices), nv))
    candidates = []
    for h, matrix in enumerate(matrices):
        fitted = [gls(matrix, response[:, v:v+1], roles, factors[v]) for v in range(nv)]
        scores[h] = [item[4][0] for item in fitted]
        candidates.append(np.column_stack([item[0][:, 0] for item in fitted]))
    hrf = np.argmin(scores, axis=0) if locked_hrf is None else np.asarray(locked_hrf)
    stage_b = np.column_stack([candidates[h][:, v] for v, h in enumerate(hrf)])
    pairs = [(i, j) for i in range(ntrial) for j in range(ntrial)
             if trial_runs[i] != trial_runs[j] and conditions[i] == conditions[j]]
    trends, stages = [], []
    for count in range(maxpc+1):
        current_roles = roles+[0]*(len(np.unique(runs))*count)
        voxel_betas = []
        for voxel, chosen_hrf in enumerate(hrf):
            matrix = augmented_design(matrices[chosen_hrf], pcs, runs, count)
            fit = gls(matrix, response[:, voxel:voxel+1], current_roles,
                      factors[voxel])
            voxel_betas.append(fit[0][:, 0])
        beta = np.column_stack(voxel_betas)
        losses = np.mean([(beta[i]-stage_b[j])**2 for i, j in pairs], axis=0)
        losses /= stage_b.var(axis=0, ddof=1)
        trends.append(np.median(losses[select_voxels]))
        stages.append(beta)
    improvement = trends[0]-np.asarray(trends)
    chosen = int(np.flatnonzero(improvement >= improvement.max()/1.05)[0])
    if fixed_count is not None: chosen = fixed_count
    final, nuisance, fitted, wins, cv = [], [], [], [], []
    for v, h in enumerate(hrf):
        matrix = augmented_design(matrices[h], pcs, runs, chosen)
        current_roles = roles+[0]*(len(np.unique(runs))*chosen)
        if len(fractions) > 1:
            loss, winner, _ = repeat_oracle(np.linalg.solve(factors[v], matrix),
                np.linalg.solve(factors[v], response[:, v:v+1]), current_roles,
                runs, trial_runs, conditions, fractions)
            index = int(winner[0])
            cv.append(loss[index, 0])
        else:
            index = 0
        beta, gamma, prediction, _, _ = gls(matrix, response[:, v:v+1],
                                          current_roles, factors[v], fractions[index])
        final.append(beta[:, 0])
        nuisance.append(gamma[:, 0])
        fitted.append(prediction[:, 0])
        wins.append(index)
    return dict(hrf=hrf, scores=scores, stage_b=stage_b, stage_c=stages[chosen],
                chosen=chosen, trend=trends, beta=np.column_stack(final),
                nuisance=np.column_stack(nuisance), fitted=np.column_stack(fitted),
                wins=wins, cv=cv)


def check_pipeline(prefix, expected, pcs, ntrial, nv, tuned=True):
    """Check intermediate choices and fits, not just the final beta endpoint."""
    data, info, labels = bucket(prefix)
    assert_allclose(data[labels['Picard_DesignIndex'], :nv], expected['hrf'])
    for h, identifier in enumerate([11, 4]):
        assert_allclose(data[labels[f'Picard_DesignMSE_ID{identifier}'], :nv], expected['scores'][h], rtol=3e-5, atol=3e-6)
    for stage in ['b', 'c']:
        first = labels[f'Stage{stage.upper()}_identity0_Coef']
        assert_allclose(data[first:first+ntrial, :nv], expected[f'stage_{stage}'], rtol=3e-5, atol=3e-5)
    assert info['PICARD_PC_COUNT'] == expected['chosen']
    assert_allclose(data[labels['Picard_PCCount'], :nv], expected['chosen'])
    actual_pcs = np.fromstring(info['PICARD_NOISE_PCS'], sep=' ').reshape(pcs.shape)
    assert_allclose(actual_pcs, pcs, atol=3e-8)
    if tuned:
        trend = [float(item.split()[1]) for item in info['PICARD_PC_TREND'].split(';')]
        assert_allclose(trend, expected['trend'], atol=2e-6)
        assert_allclose(data[labels['Picard_CV_Index'], :nv], expected['wins'])
        assert_allclose(data[labels['Picard_CV_MSE'], :nv], expected['cv'], rtol=4e-5, atol=4e-5)
    assert_allclose(data[:ntrial, :nv], expected['beta'], atol=4e-5, rtol=4e-5)
    assert_allclose(read(prefix, 'nuisance')[:, :nv], expected['nuisance'], atol=4e-5, rtol=4e-5)
    assert_allclose(read(prefix, 'fitted')[:, :nv], expected['fitted'], atol=4e-5, rtol=4e-5)
    assert info['PICARD_DESIGN_DIAGNOSTIC_START'] == labels['Picard_DesignID']
    assert 'independent_outer_evaluation=0' in info['PICARD_STAGE_CONFIG']
    assert info['PICARD_HRF_CRITERION'] == 'weighted_residual_SSE_over_n'
    assert np.all(data[:ntrial, nv:] == 0)
    assert np.all(data[labels['Picard_PCCount'], nv:] == -1)
    return data, info, labels


def main():
    with tempfile.TemporaryDirectory(prefix='picard-ng2-') as work:
        os.chdir(work)
        rng = np.random.default_rng(709)
        nr, nt, nv, selected, ntrial = 3, 36, 14, 6, 9
        nfull = nr*nt
        run_full = np.repeat(np.arange(nr), nt)
        good = np.delete(np.arange(nfull), [3, 4, 40, 76])
        runs = run_full[good]
        trial_runs = np.repeat(np.arange(nr), 3)
        conditions = np.tile(np.arange(3), nr)
        roles = [1]*ntrial+[0]*nr
        designs = []
        for h in range(2):
            matrix = np.zeros((nfull, ntrial+nr))
            for r in range(nr):
                matrix[run_full == r, r*3:r*3+3] = rng.normal(size=(nt, 3))
                matrix[run_full == r, ntrial+r] = 1
            designs.append(matrix)
        noise = rng.normal(size=(nfull, 2))
        for r in range(nr): noise[run_full == r, 0] += (r-1)*3*designs[0][run_full == r, r*3]
        response = noise@rng.normal(size=(2, nv))*2
        truth = np.tile([1, 3, -2], nr)*3
        for v in range(selected): response[:, v] += designs[v % 2][:, :ntrial]@truth
        response += rng.normal(size=response.shape)*.1
        response = response.astype(np.float32).astype(float)
        manifest, matrices = candidate_manifest('hrf', [x[good] for x in designs], roles,
            good, np.arange(nr)*nt, ids=(11, 4), trial_runs=trial_runs, conditions=conditions)
        volume('data.nii', response)
        volume('pool.nii', [np.arange(nv) >= selected])
        volume('analysis.nii', [np.arange(nv) < selected])
        volume('select.nii', [np.arange(nv) < 4])
        np.savetxt('trials.1D', np.c_[np.arange(ntrial), np.arange(ntrial), conditions, trial_runs], fmt='%d')
        pilot = xmat('pilot.1D', matrices[0], roles, good, nfull, np.arange(nr)*nt)
        pcs = pca_reference(response[good], runs, np.arange(selected, nv), 2)
        base = ['-input', 'data.nii', '-candidate_manifest', manifest, '-candidate_policy', 'hrf',
                '-trial_info', 'trials.1D', '-single_trial', '-noise_pool', 'pool.nii',
                '-mask', 'analysis.nii', '-pc_counts', '0,1,2', '-pc_mask', 'select.nii',
                '-fractions', '1,.7,.3', '-save_stages', '-save_candidates', '-save_nuisance',
                '-save_fitted', '-save_residuals']
        # Distinct pairs include a group containing no PC-selection voxels.
        # A naive group-local median would either fail or pick a different count.
        parameters = np.array([[.6,.2],[.8,-.5],[-.4,.2],[.6,.2],[.8,-.5],[.3,-.2]], np.float32).astype(float)
        all_parameters = np.zeros((nv, 2)); all_parameters[:selected] = parameters
        volume('cov.nii', all_parameters.T)
        for gap in ['preserve', 'reset']:
            expected = staged_reference(matrices, response[good, :selected], roles, good, runs,
                parameters, 'arma11', gap, pcs, trial_runs, conditions, [1,.7,.3], np.arange(4))
            assert expected['chosen'] > 0
            prefix = 'mapped_'+gap
            opts = ['-covariance', 'arma11', '-covariance_map', 'cov.nii', '-covariance_gaps', gap]
            run([BINARY, *base, *opts, '-batch_size', '2', '-prefix', prefix])
            data, info, labels = check_pipeline(prefix, expected, pcs, ntrial, selected)
            assert 'groups=4 ' in info['PICARD_COVARIANCE_GROUPS']
            run([BINARY, *base, *opts, '-batch_size', '5', '-jobs',
                 os.environ.get('PICARD_TEST_JOBS', '2'), '-prefix', 'parallel_'+gap])
            assert_allclose(bucket('parallel_'+gap)[0], data, atol=3e-5)
            assert_allclose(read(prefix, 'residuals')[:, :selected]+read(prefix, 'fitted')[:, :selected],
                            response[good, :selected], atol=4e-5)
        print('PASS shared PC selection across mapped covariance groups, both gaps and scheduling')
        # Native pilot estimates are made once and frozen before HRF selection.
        # The independent QR estimator has no access to Picard's reported maps.
        native, _ = reml_grid(pilot, response[good, :selected], good, runs, 'arma11')
        expected = staged_reference(matrices, response[good, :selected], roles, good, runs,
            native[:, :2], 'arma11', 'preserve', pcs, trial_runs, conditions, [1,.7,.3], np.arange(4))
        native_opts = ['-covariance', 'arma11', '-covariance_estimate', 'pilot.1D', '-covariance_freeze']
        run([BINARY, *base, *native_opts, '-prefix', 'native'])
        data, info, labels = check_pipeline('native', expected, pcs, ntrial, selected)
        assert_allclose(data[-2:, :selected].T, native[:, :2], atol=1e-6)
        assert 'native_freeze_opt_in=1' in info['PICARD_COVARIANCE']
        all_parameters[:selected] = native[:, :2]
        volume('native_params.nii', all_parameters.T)
        run([BINARY, *base, '-covariance', 'arma11', '-covariance_map', 'native_params.nii', '-prefix', 'imported'])
        imported, _, ilabels = bucket('imported')
        for name, at in ilabels.items():
            if name: assert_allclose(imported[at], data[labels[name]], atol=4e-5, rtol=4e-5)
        # Explicit white covariance is an ablation of all temporal weighting.
        run([BINARY, *base, '-covariance', 'white', '-prefix', 'white'])
        run([BINARY, *base, '-prefix', 'ordinary'])
        white, _, wl = bucket('white'); ordinary, _, ol = bucket('ordinary')
        for name, at in ol.items():
            if name: assert_allclose(ordinary[at], white[wl[name]], atol=3e-5)
        print('PASS native frozen pilot against QR oracle, map reuse and white ablation')
        # Fixing count/penalty ablates both selection steps. Zero PCs must
        # coincide with ordinary HRF-candidate GLS; at fraction 1, final betas
        # must coincide with the selected-PC unpenalized intermediate stage.
        fixed = base.copy()
        at = fixed.index('-pc_counts'); fixed[at:at+4] = ['-noise_pcs', '0']
        at = fixed.index('-fractions'); fixed[at:at+2] = ['-fraction', '.7']
        for count, fraction in [(0, .7), (1, 1.)]:
            fixed[fixed.index('-noise_pcs')+1] = str(count)
            fixed[fixed.index('-fraction')+1] = str(fraction)
            run([BINARY, *fixed, *native_opts, '-prefix', f'fixed{count}'])
            current = staged_reference(matrices, response[good, :selected], roles, good, runs,
                native[:, :2], 'arma11', 'preserve', pcs, trial_runs, conditions,
                [fraction], np.arange(4), fixed_count=count)
            values, header, names = bucket(f'fixed{count}')
            assert_allclose(values[:ntrial, :selected], current['beta'], atol=4e-5, rtol=4e-5)
            assert header['PICARD_PC_COUNT'] == count
            if fraction == 1:
                first = names['StageC_identity0_Coef']
                assert_allclose(values[:ntrial], values[first:first+ntrial], atol=2e-6)
        ordinary_args = ['-input','data.nii','-candidate_manifest',manifest,
                         '-candidate_policy','hrf','-fraction','.7','-mask','analysis.nii']
        run([BINARY, *ordinary_args, *native_opts, '-prefix','no_pc_ordinary'])
        assert_allclose(read('no_pc_ordinary','betas'), read('fixed0','betas'), atol=4e-5)
        # Prediction mode reuses the same frozen covariance and ordinary fold
        # machinery. This deliberately tests conditional prediction; it does
        # not claim the full-data pilot was trained independently of the folds.
        xmat('prediction.1D', matrices[0], roles, good, nfull, np.arange(nr)*nt)
        predict = ['-input','data.nii','-matrix','prediction.1D','-mask','analysis.nii',
                   '-cv_predict','-cv_block_size','12','-cv_nuisance','project',
                   '-fractions','1,.7,.3','-save_cv']
        run([BINARY, *predict, *native_opts, '-prefix','prediction_native'])
        run([BINARY, *predict, '-covariance','arma11','-covariance_map','native_params.nii',
             '-prefix','prediction_imported'])
        predicted, _, plabels = bucket('prediction_native')
        imported, _, ilabels = bucket('prediction_imported')
        for name, at in ilabels.items():
            if name: assert_allclose(imported[at], predicted[plabels[name]], atol=5e-5, rtol=5e-5)
        # The covariance grouping path also supports ordinary prediction-based
        # candidate selection, which is separate from the staged HRF policy.
        candidate_predict = predict.copy()
        at = candidate_predict.index('-matrix')
        candidate_predict[at:at+2] = ['-candidate_manifest', manifest,
                                     '-candidate_policy', 'predict']
        run([BINARY, *candidate_predict, *native_opts, '-prefix','predict_candidates_native'])
        run([BINARY, *candidate_predict, '-covariance','arma11',
             '-covariance_map','native_params.nii','-prefix','predict_candidates_imported'])
        assert_allclose(read('predict_candidates_native','betas'),
                        read('predict_candidates_imported','betas'), atol=5e-5, rtol=5e-5)
        # The coarse pool criterion is explicitly unweighted original-time R2,
        # even when subsequent trial fits use voxel-specific covariance.
        auto_pool = base.copy()
        at = auto_pool.index('-noise_pool')
        auto_pool[at:at+2] = ['-noise_model','pilot.1D','-noise_search','pool.nii',
                             '-noise_r2','1']
        run([BINARY, *auto_pool, *native_opts, '-prefix','auto_pool'])
        automatic, automatic_info, auto_labels = bucket('auto_pool')
        assert_allclose(automatic[:ntrial], data[:ntrial], atol=4e-5)
        assert automatic_info['PICARD_PC_COUNT'] == info['PICARD_PC_COUNT']
        assert_allclose(np.fromstring(automatic_info['PICARD_NOISE_PCS'], sep=' '),
                        pcs.ravel(), atol=3e-8)
        # Appended PCs stay unpenalized for the other existing solvers too.
        # Independent weighted augmented least squares checks banded ridge;
        # the existing sklearn oracle checks the elastic-net objective.
        np.savetxt('groups.1D', [list(np.arange(ntrial)%2)+[-1]*nr], fmt='%d')
        np.savetxt('bands.1D', [[.1,.8]])
        penalized = fixed.copy()
        at = penalized.index('-fraction'); del penalized[at:at+2]
        for kind, options, parameter, groups in [
                ('banded', ['-penalty','banded','-band_groups','groups.1D',
                            '-band_penalties','bands.1D'], [.1,.8], np.arange(ntrial)%2),
                ('enet', ['-penalty','enet','-lambda','.15','-l1_ratio','.7'], .15, None)]:
            run([BINARY, *penalized, *native_opts, *options, '-prefix',kind])
            actual = read(kind, 'betas')
            for voxel, hrf in enumerate(expected['hrf']):
                matrix = augmented_design(matrices[hrf], pcs, runs, 1)
                lower = np.linalg.cholesky(correlation(good, runs, 'arma11', *native[voxel,:2]))
                reference = phase6_fit(np.linalg.solve(lower,matrix),
                    np.linalg.solve(lower,response[good,voxel:voxel+1]),
                    roles+[0]*nr, good, runs, 0, kind, parameter, groups)
                assert_allclose(actual[:,voxel], reference[0][:,0], atol=6e-5, rtol=6e-5)
        print('PASS fixed-count, zero-PC, unpenalized and frozen-pilot prediction ablations')
        # A simple adversarial HRF example distinguishes weighted from raw MSE.
        # Slow residuals are cheaper under positive AR correlation. Consequently
        # GLS prefers the fast predictor even though raw error prefers the slow one.
        n = 32; t = np.arange(n)
        slow, fast = np.cos(2*np.pi*t/n), (-1.)**t/np.sqrt(2)
        hrf_manifest, choices = candidate_manifest('scoring',
            [np.c_[slow, np.ones(n)], np.c_[fast, np.ones(n)]], [1,0], t, [0], ids=(11,4))
        sample = (1.2*slow+fast)[:, None]
        volume('scoring.nii', sample)
        lower = np.linalg.cholesky(correlation(t, np.zeros(n,int), 'ar1', .8))
        fits = [gls(x, sample, [1,0], lower) for x in choices]
        weighted_winner = int(np.argmin([fit[4][0] for fit in fits]))
        raw_winner = int(np.argmin([np.mean(fit[3]**2) for fit in fits]))
        assert weighted_winner != raw_winner
        run([BINARY, '-input','scoring.nii','-candidate_manifest',hrf_manifest,
             '-candidate_policy','hrf','-fraction','1','-covariance','ar1',
             '-covariance_ar','.8','-prefix','scoring'])
        actual, _, names = bucket('scoring')
        assert actual[names['Picard_DesignIndex'], 0] == weighted_winner
        assert_allclose(actual[names['Picard_DesignMSE'], 0], fits[weighted_winner][4][0], atol=2e-6)
        # Adaptive native use is always explicit; failure leaves no output.
        result = run([BINARY, *base, *native_opts[:-1], '-prefix', 'rejected'], ok=False)
        assert 'covariance_freeze' in result
        assert not list(Path('.').glob('rejected*'))
        print('PASS covariance-weighted HRF criterion and explicit adaptive opt-in')


if __name__ == '__main__':
    main()
