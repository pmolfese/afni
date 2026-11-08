#!/usr/bin/env python3
"""Staged single-trial integration: independent NumPy fits/PCA and failures.

Reuse only fixture I/O and the independent least-squares oracle from the
pre-existing integration suite. No network or GLMsingle installation needed.
"""
from pathlib import Path
source = Path(__file__).with_name('test_picard.py').read_text()
exec(compile(source.split("with tempfile.TemporaryDirectory(prefix='picard-tests-')")[0],
             str(Path(__file__).with_name('test_picard.py')), 'exec'))


def bucket(prefix):
    img = nib.load(prefix+'+orig.HEAD')
    return np.asarray(img.dataobj).reshape(img.shape[0], -1).T, img.header.info


def label_map(info):
    return {label: i for i, label in enumerate(info['BRICK_LABS'].split('~'))}


with tempfile.TemporaryDirectory(prefix='picard-single-') as folder:
    os.chdir(folder)
    rng = np.random.default_rng(709)
    nr, nt, nv, p = 3, 36, 14, 9
    n = nr*nt
    rows = np.repeat(np.arange(nr), nt)
    trials = np.repeat(np.arange(nr), 3)
    conditions = np.tile(np.arange(3), nr)
    roles = [1]*p+[0]*nr
    matrices = []
    for h in range(2):
        x = np.zeros((n, p+nr))
        for r in range(nr):
            x[rows == r, r*3:r*3+3] = rng.normal(size=(nt, 3))
            x[rows == r, p+r] = 1
        matrices.append(x)
    manifest, matrices = candidate_manifest('hrf', matrices, roles, np.arange(n),
                  np.arange(nr)*nt, ids=(11, 4), trial_runs=trials,
                  conditions=conditions)
    # Strong shared temporal noise with unequal run-dependent task correlation.
    noise = rng.normal(size=(n, 2))
    for r in range(nr):
        noise[rows == r, 0] += (r-1)*3*matrices[0][rows == r, r*3]
    y = noise @ rng.normal(size=(2, nv))*2
    winners = np.arange(nv) % 2
    truth = np.tile(np.array([1, 3, -2]), nr)
    for v in range(6):
        y[:, v] += matrices[winners[v]][:, :p]@truth*3
    y += rng.normal(size=y.shape)*.1
    y = y.astype(np.float32).astype(float)
    volume('data.nii', y)
    volume('pool.nii', [np.arange(nv) >= 6])
    volume('analysis.nii', [np.arange(nv) < 6])
    volume('select.nii', [np.arange(nv) < 4])
    np.savetxt('trials.1D', np.c_[np.arange(p), np.arange(p), conditions, trials], fmt='%d')
    args = ['-input','data.nii','-candidate_manifest',manifest,
            '-candidate_policy','hrf','-trial_info','trials.1D','-single_trial',
            '-noise_pool','pool.nii','-mask','analysis.nii']
    selection = ['-pc_counts','0,1,2','-pc_mask','select.nii']
    run([BINARY,*args,*selection,'-fractions','1,.7,.3','-save_stages',
         '-save_nuisance','-save_fitted','-save_residuals','-save_candidates',
         '-save_cv','-prefix','staged','-batch_size','2'])
    data, info = bucket('staged'); labels = label_map(info)
    # PCA is deliberately independent: SVD of normalized pool responses,
    # rather than SVD of their temporal Gram matrix as in production.
    pcs = np.zeros((n, 2))
    for r in range(nr):
        pool = y[rows == r, 6:].copy()
        pool -= pool.mean(axis=0)
        pool /= np.linalg.norm(pool, axis=0)
        u, _, _ = np.linalg.svd(pool, full_matrices=False)
        u = u[:, :2]/u[:, :2].std(axis=0, ddof=1)
        u *= np.sign(u[np.argmax(np.abs(u), axis=0), np.arange(2)])
        pcs[rows == r] = u
    actual_pcs = np.fromstring(info['PICARD_NOISE_PCS'], sep=' ').reshape(n, 2)
    assert_allclose(actual_pcs, pcs, atol=2e-8)
    fit_b = [oracle(x, y[:, :6], roles, 1) for x in matrices]
    hrf = np.argmin([np.mean(f[3]**2, axis=0) for f in fit_b], axis=0)
    assert set(hrf) == {0, 1}
    assert_allclose(data[labels['Picard_DesignIndex'], :6], hrf)
    reference = np.column_stack([fit_b[h][0][:, v] for v, h in enumerate(hrf)])
    assert_allclose(data[labels['StageB_identity0_Coef']:][:p, :6], reference, atol=2e-5)
    augmented, losses, fits = [], [], []
    pairs = [(i,j) for i in range(p) for j in range(p)
             if trials[i] != trials[j] and conditions[i] == conditions[j]]
    for count in range(3):
        extra = np.zeros((n, nr*count))
        for r in range(nr): extra[rows == r, r*count:(r+1)*count] = pcs[rows == r, :count]
        design = [np.c_[x, extra] for x in matrices]
        current = [oracle(x, y[:, :6], roles+[0]*(nr*count), 1)[0] for x in design]
        beta = np.column_stack([current[h][:, v] for v,h in enumerate(hrf)])
        loss = np.mean([(beta[i]-reference[j])**2 for i,j in pairs], axis=0)/reference.var(axis=0, ddof=1)
        losses.append(np.median(loss[:4])); augmented.append(design); fits.append(beta)
    improvement = losses[0]-np.asarray(losses)
    chosen = int(np.flatnonzero(improvement >= improvement.max()/1.05)[0])
    # Fixture exercises a nonzero selected count, not just the no-PC endpoint.
    assert chosen > 0, losses
    assert info['PICARD_PC_COUNT'] == chosen
    trend = [float(line.split()[1]) for line in info['PICARD_PC_TREND'].split(';')]
    assert_allclose(trend, losses, atol=2e-7)
    assert_allclose(data[labels['StageC_identity0_Coef']:][:p, :6], fits[chosen], atol=2e-5)
    chosen_roles = roles+[0]*(nr*chosen)
    for h in range(2):
        vox = np.flatnonzero(hrf == h)
        x = augmented[chosen][h]
        loss, win, _ = repeat_oracle(x, y[:, vox], chosen_roles, rows, trials, conditions, [1,.7,.3])
        assert_allclose(data[labels['Picard_CV_Index'], vox], win)
        for o,v in enumerate(vox):
            expected = oracle(x,y[:, v:v+1],chosen_roles,[1,.7,.3][win[o]])
            assert_allclose(data[:p, v], expected[0][:,0], atol=3e-5)
            assert_allclose(read('staged','fitted')[:,v], expected[2][:,0], atol=3e-5)
    assert info['PICARD_DIMENSIONS'][3] == nr*(1+chosen)
    assert len(info['PICARD_COLUMN_ROLES']) == p+nr*(1+chosen)
    assert 'NOT independent test performance' in info['PICARD_SINGLE_TRIAL']
    assert np.all(data[labels['Picard_PCCount'], 6:] == -1)
    assert np.all(data[:p, 6:] == 0)
    print('PASS staged HRF/PCA/PC selection/final repeat fits and bucket provenance')
    # Batch/thread changes preserve upstream selection and final fits.
    run([BINARY,*args,*selection,'-fractions','1,.7,.3','-save_stages',
         '-save_nuisance','-save_fitted','-save_residuals','-save_candidates',
         '-save_cv','-prefix','parallel','-batch_size','5','-jobs',os.environ.get('PICARD_TEST_JOBS','2')])
    assert_allclose(bucket('parallel')[0], data, atol=2e-5)
    # Fixed counts support one supplied design, and sparse/support/AR1 adapters.
    fixed = ['-input','data.nii','-matrix',str(Path(manifest).parent/'design 11.xmat.1D'),
             '-single_trial','-trial_info','trials.1D','-noise_pool','pool.nii',
             '-mask','analysis.nii','-noise_pcs','1']
    run([BINARY,*fixed,'-fraction','.7','-prefix','fixed','-save_stages'])
    expected = oracle(augmented[1][0], y[:,:6], roles+[0]*nr, .7)[0]
    assert_allclose(read('fixed','betas')[:,:6], expected, atol=2e-5)
    run([BINARY,*fixed,'-penalty','enet','-lambda','.1','-l1_ratio','.5',
         '-support_refit','-prewhite_ar1','.2','-prefix','sparse'])
    assert np.isfinite(bucket('sparse')[0]).all()
    # Explicit model-R2 pool must reproduce an explicit mask when its threshold
    # includes all search voxels. The coarse design shares sampling/nuisance.
    model = xmat('coarse.xmat.1D',matrices[0],[1]*p+[0]*nr,runs=np.arange(nr)*nt)
    auto = args.copy(); at=auto.index('-noise_pool'); auto[at:at+2]=[
        '-noise_model','coarse.xmat.1D','-noise_search','pool.nii','-noise_r2','1']
    run([BINARY,*auto,'-noise_pcs','1','-fraction','.7','-prefix','auto'])
    auto_data, auto_info = bucket('auto'); auto_labels=label_map(auto_info)
    fit = oracle(model,y[:,6:],roles,1)
    nuisance = model[:,p:]; residual=y[:,6:]-nuisance@np.linalg.lstsq(nuisance,y[:,6:],rcond=None)[0]
    r2=1-np.sum(fit[3]**2,axis=0)/np.sum(residual**2,axis=0)
    assert_allclose(auto_data[auto_labels['Picard_NoiseModelR2'],6:],r2,atol=2e-6)
    assert_allclose(np.fromstring(auto_info['PICARD_NOISE_PCS'],sep=' ').reshape(n,1)[:,0],pcs[:,0],atol=2e-8)
    print('PASS batching, fixed count, sparse/AR1 adapter, and coarse-model pool')
    good=np.delete(np.arange(n),[3,4,40,76])
    cruns=rows[good]
    cx=xmat('censored.xmat.1D',matrices[0][good],roles,good,n,np.arange(nr)*nt)
    cp=np.zeros((len(good),nr))
    for r in range(nr):
        local=y[good][:,6:][cruns==r].copy()
        local-=local.mean(axis=0);local/=np.linalg.norm(local,axis=0)
        u=np.linalg.svd(local,full_matrices=False)[0][:,0]
        cp[cruns==r,r]=u/u.std(ddof=1)
    cx=np.c_[cx,cp]
    censored=fixed.copy();censored[censored.index(str(Path(manifest).parent/'design 11.xmat.1D'))]='censored.xmat.1D'
    run([BINARY,*censored,'-fraction','.7','-prewhite_ar1','.4',
         '-save_fitted','-save_residuals','-prefix','censored'])
    expected=phase6_fit(cx,y[good,:6],np.array(roles+[0]*nr),good,cruns,.4,'ridge',.7,None)
    assert_allclose(read('censored','betas')[:,:6],expected[0],atol=3e-5)
    assert_allclose(read('censored','fitted')[:,:6],expected[2],atol=3e-5)
    # Unit-marginal AR1 with the same reset policy differs by an overall
    # whitening scale: fractional ridge must preserve betas after PC addition.
    run([BINARY,*censored,'-fraction','.7','-covariance','ar1',
         '-covariance_ar','.4','-covariance_gaps','reset',
         '-save_fitted','-prefix','covariance_staged'])
    assert_allclose(read('covariance_staged','betas')[:,:6],expected[0],atol=3e-5)
    assert_allclose(read('covariance_staged','fitted')[:,:6],expected[2],atol=3e-5)
    # Banded ridge keeps its group identities when PCs are appended to Z.
    np.savetxt('groups.1D',[[*(np.arange(p)%2),*([-1]*nr)]],fmt='%d')
    np.savetxt('bands.1D',[[.1,.8]])
    run([BINARY,*fixed,'-penalty','banded','-band_groups','groups.1D',
         '-band_penalties','bands.1D','-prefix','bands'])
    expected=phase6_fit(augmented[1][0],y[:,:6],np.array(roles+[0]*nr),
                       np.arange(n),rows,0,'banded',np.array([.1,.8]),np.arange(p)%2)
    assert_allclose(read('bands','betas')[:,:6],expected[0],atol=3e-5)
    print('PASS censored-row PCA/AR1 and appended-PC banded-ridge identities')
    # Failures must leave no partial output datasets.
    def reject(extra, message, base=args):
        result=run([BINARY,*base,*extra,'-prefix','failure'],ok=False)
        assert message in result, result
        assert not list(Path('.').glob('failure*'))
    reject(['-noise_pcs','35','-fraction','1'],'rank')
    reject(['-pc_counts','1,2','-pc_mask','select.nii','-fraction','1'],'zero')
    reject(['-pc_counts','0,1','-pc_mask','pool.nii','-fraction','1'],'pc_mask')
    reject(['-noise_pcs','1','-fraction','1','-memory_mb','.001'],'memory')
    volume('empty.nii',[np.zeros(nv)])
    bad=args.copy(); bad[bad.index('pool.nii')]='empty.nii'
    reject(['-noise_pcs','1','-fraction','1'],'empty',bad)
    run([BINARY,*args,'-noise_pcs','0','-fraction','1','-prefix','zero'])
    assert_allclose(read('zero','betas')[:,:6],reference,atol=2e-5)
    zero_y=y.copy(); zero_y[:,:6]=0
    volume('undefined.nii',zero_y)
    bad=args.copy(); bad[bad.index('data.nii')]='undefined.nii'
    reject([*selection,'-fraction','1'],'no defined',bad)
    bad=args.copy(); at=bad.index('-noise_pool'); bad[at:at+2]=[
        '-noise_model','coarse.xmat.1D','-noise_search','pool.nii','-noise_r2','0']
    reject(['-noise_pcs','1','-fraction','1'],'empty',bad)
    # A failing companion write removes every dataset from this transaction.
    Path('blocked').write_text('preserve')
    reject(['-noise_pcs','1','-fraction','1','-fitts','blocked/out'], 'cannot')
    assert Path('blocked').read_text()=='preserve'
    print('PASS staged failure contracts and zero-PC endpoint')
    ref=json.loads(Path(__file__).with_name('picard_single_reference.json').read_text())
    rx=np.asarray(ref['matrix']);ry=np.asarray(ref['data']);rp=len(ref['conditions'])
    rn=len(ry);rnt=ref['run_length'];rnr=rn//rnt
    xmat('pinned.xmat.1D',rx,[1]*rp+[0]*rnr,runs=np.arange(rnr)*rnt)
    volume('pinned.nii',ry)
    volume('pinned_pool.nii',[np.arange(ry.shape[1])>=3])
    volume('pinned_mask.nii',[np.arange(ry.shape[1])<3])
    np.savetxt('pinned_trials.1D',np.c_[np.arange(rp),np.arange(rp),
               ref['conditions'],ref['trial_runs']],fmt='%d')
    run([BINARY,'-input','pinned.nii','-matrix','pinned.xmat.1D',
         '-trial_info','pinned_trials.1D','-single_trial',
         '-noise_pool','pinned_pool.nii','-mask','pinned_mask.nii',
         '-pc_mask','pinned_mask.nii','-pc_counts','0,1,2',
         '-fraction','1','-save_stages','-save_fitted','-prefix','pinned'])
    out,header=bucket('pinned');names=label_map(header)
    assert header['PICARD_PC_COUNT']==ref['chosen']
    assert_allclose(np.fromstring(header['PICARD_NOISE_PCS'],sep=' ').reshape(rn,2),
                    ref['pcs'],atol=3e-6)
    expected=np.median(np.asarray(ref['badness']),axis=0)/ref['pair_count']
    actual=[float(line.split()[1]) for line in header['PICARD_PC_TREND'].split(';')]
    assert_allclose(actual,expected,atol=3e-6)
    assert_allclose(out[:rp,:3],ref['final_beta'],atol=3e-5)
    assert_allclose(out[names['StageB_col0_Coef']:][:rp,:3],ref['stage_b'],atol=3e-5)
    assert_allclose(read('pinned','fitted')[:,:3],ref['final_fitted'],atol=3e-5)
    print('PASS pinned GLMsingle PCA/scoring/stop and controlled staged fixture:',ref['commit'])
print('PASS 3dPicard staged single-trial integration')
