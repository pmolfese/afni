#!/usr/bin/env python3
"""Picard fixed/tuned fits: independent numerical oracles and CLI contracts.

Requires numpy, scipy, nibabel and scikit-learn. Existing 3dDeconvolve is a fixture producer,
never a production dependency of 3dPicard. Fracridge fixture provenance is
recorded in picard_fracridge_reference.json. The repeated-condition scoring
fixture pins GLMsingle independently; no network is used by these tests.
"""
import csv
import gzip
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import nibabel as nib
import numpy as np
from numpy.testing import assert_allclose
from scipy.optimize import brentq

BINARY = str(Path(sys.argv[1]).resolve())
DECONVOLVE = str(Path(sys.argv[2]).resolve())
FIXTURE = Path(__file__).with_name('picard_fracridge_reference.json')
ENV = dict(os.environ, AFNI_NIFTI_TYPE_WARN='NO', AFNI_NIFTI_VIEW='orig',
           OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', OMP_NUM_THREADS='2',
           AFNI_WRITE_1D_AS_PREFIX='YES', AFNI_FLOATSCAN='YES')

def run(args, ok=True):
    p = subprocess.run([str(a) for a in args], env=ENV, text=True, capture_output=True)
    if (p.returncode == 0) != ok:
        raise AssertionError(f'{args}\nreturn={p.returncode}\n{p.stdout}\n{p.stderr}')
    return p.stdout + p.stderr

def volume(name, values, dtype=np.float32, slope=None):
    values = np.asarray(values, dtype=dtype)
    # rows=time, columns=voxel; simple volume with one spatial axis.
    img = nib.Nifti1Image(values.T[:, None, None, :], np.eye(4))
    if slope is not None:
        img.header.set_slope_inter(slope, 0)
    nib.save(img, name)

def read(prefix, kind):
    suffix={'fitted':'_fitts','residuals':'_errts'}.get(kind,'')
    img=nib.load(f'{prefix}{suffix}+orig.HEAD')
    data=np.asarray(img.dataobj).reshape(img.shape[0],-1).T
    if suffix:return data
    info=img.header.info
    p=info['PICARD_DIMENSIONS'][2];diag=info['PICARD_DIAGNOSTIC_START']
    if kind=='betas':return data[:p]
    if kind=='nuisance':return data[p:diag]
    if kind=='diagnostics':return data[diag:]
    raise ValueError(kind)

def xmat(name, x, groups=None, good=None, nfull=None, runs=None, labels=None):
    x = np.asarray(x, dtype=np.float32)
    attrs = {'ColumnLabels': ' ; '.join(labels or [f'col{j}' for j in range(x.shape[1])])}
    if groups is not None: attrs['ColumnGroups'] = ','.join(map(str, groups))
    if good is not None: attrs['GoodList'] = ','.join(map(str, good))
    if nfull is not None: attrs['NRowFull'] = str(nfull)
    if runs is not None: attrs['RunStart'] = ','.join(map(str, runs))
    with open(name, 'w') as f:
        f.write(f'# <matrix\n# ni_type = "{x.shape[1]}*float"\n# ni_dimen = "{len(x)}"\n')
        for key, val in attrs.items(): f.write(f'# {key} = "{val}"\n')
        f.write('# >\n')
        np.savetxt(f, x, fmt='%.9g')
        f.write('# </matrix>\n')
    return x.astype(float)

def oracle(matrix, y, roles, fraction):
    z, x = matrix[:, np.asarray(roles) == 0], matrix[:, np.asarray(roles) == 1]
    pz = np.linalg.pinv(z, rcond=1e-10) if z.shape[1] else np.empty((0, len(z)))
    xp = x-z@(pz@x); yp = y-z@(pz@y)
    for a, orig in [(xp,x), (yp,y)]:
        tiny = np.linalg.norm(a, axis=0) <= 64*np.finfo(float).eps*max(z.shape)*np.linalg.norm(orig, axis=0)
        a[:,tiny] = 0
    reference = np.linalg.lstsq(xp, yp, rcond=1e-10)[0]
    beta = np.zeros_like(reference); penalties = []; norms = np.linalg.norm(reference, axis=0)
    for o, norm in enumerate(norms):
        if norm == 0: penalty=0
        elif fraction == 0: penalty=np.inf
        elif fraction == 1: penalty=0; beta[:,o]=reference[:,o]
        else:
            # Independent augmented least-squares + Brent root (not production's
            # rotated-coefficient norm and bounded bisection implementation).
            def solve(pen):
                return np.linalg.lstsq(np.vstack([xp, np.sqrt(pen)*np.eye(x.shape[1])]),
                                      np.r_[yp[:,o], np.zeros(x.shape[1])], rcond=1e-14)[0]
            hi=np.linalg.norm(xp, 2)**2*(1/fraction-1)*1.01
            penalty=brentq(lambda v: np.linalg.norm(solve(v))/norm-fraction, 0, hi,
                          xtol=1e-13, rtol=1e-13)
            beta[:,o]=solve(penalty)
        penalties.append(penalty)
    gamma=pz@(y-x@beta)
    fitted=x@beta+z@gamma
    return beta,gamma,fitted,y-fitted,np.array(penalties),norms

def check_fit(prefix, matrix, y, roles, fraction, mask=None):
    b,g,f,r,pen,norm=oracle(matrix,y,roles,fraction)
    keep=np.ones(y.shape[1],bool) if mask is None else np.asarray(mask,bool)
    assert_allclose(read(prefix,'betas')[:,keep], b[:,keep], rtol=3e-5, atol=3e-6)
    if np.any(np.asarray(roles)==0):
        assert_allclose(read(prefix,'nuisance')[:,keep], g[:,keep], rtol=3e-5, atol=3e-6)
    assert_allclose(read(prefix,'fitted')[:,keep], f[:,keep], rtol=3e-5, atol=3e-6)
    assert_allclose(read(prefix,'residuals')[:,keep], r[:,keep], rtol=3e-5, atol=3e-6)
    dg=read(prefix,'diagnostics')
    assert_allclose(dg[2,keep], np.where(np.isinf(pen),-1,pen)[keep], rtol=3e-5, atol=3e-6)
    assert_allclose(dg[1,keep], np.where(np.isinf(pen),-1,pen/len(matrix))[keep], rtol=3e-5, atol=3e-6)
    assert_allclose(dg[3,keep], norm[keep], rtol=3e-5, atol=3e-6)
    for o in np.where(keep)[0]:
        if norm[o]>0: assert_allclose(dg[4,o],fraction,atol=1e-7)
        else: assert dg[4,o]==-1 and dg[6,o]==1
    if not keep.all():
        assert np.all(read(prefix,'betas')[:,~keep]==0) and np.all(dg[6,~keep]==-1)

def failure_case(args, text=None, prefix='bad'):
    msg=run([BINARY,*args,'-prefix',prefix],ok=False)
    if text: assert text in msg, msg
    assert not (list(Path('.').glob(prefix+'.*'))+list(Path('.').glob(prefix+'+*'))), 'failure left partial outputs'

def repeat_oracle(x, y, roles, row_runs, trial_runs, conditions, fractions):
    """Direct per-fold augmented least squares, never a full-data ridge path."""
    roles = np.asarray(roles)
    penalized = np.flatnonzero(roles)
    reference = oracle(x, y, roles, 1)[0]
    loss = np.zeros((len(fractions)+1, y.shape[1]))
    pairs = 0
    for heldout in np.unique(row_runs):
        train_rows = row_runs != heldout
        train_trials = np.flatnonzero(trial_runs != heldout)
        # Drop held-out (and all-zero) columns, retaining original order.
        active = np.any(x[train_rows] != 0, axis=0)
        columns = np.flatnonzero(active)
        assert np.array_equal(columns[roles[active] == 1], penalized[train_trials])
        train_x = x[train_rows][:, active]
        for ci, fraction in enumerate([*fractions, 1]):
            train_beta = oracle(train_x, y[train_rows], roles[active], fraction)[0]
            for test in np.flatnonzero(trial_runs == heldout):
                matching = np.flatnonzero(conditions[train_trials] == conditions[test])
                differences = train_beta[matching]-reference[test]
                loss[ci] += np.sum(differences**2, axis=0)
                if ci == 0: pairs += len(matching)
    loss /= pairs
    # Larger fractions win exact ties, regardless of input grid ordering.
    ordered = sorted(range(len(fractions)), key=lambda i: -fractions[i])
    winners = np.array([min(ordered, key=lambda i: loss[i, v])
                        for v in range(y.shape[1])])
    return loss, winners, pairs


def repeat_case(tag, matrix, y, roles, row_runs, trial_runs, conditions,
                fractions, good=None, nfull=None, starts=None, existing=None):
    """Check scores, chosen fractions, final fits, identities and scheduling."""
    matrix = np.asarray(matrix, np.float32).astype(float)
    y = np.asarray(y, np.float32).astype(float)
    if good is None: good = np.arange(len(matrix))
    if nfull is None: nfull = len(matrix)
    if starts is None: starts = np.flatnonzero(np.r_[True, np.diff(row_runs)!=0])
    filename = existing or f'{tag}.xmat.1D'
    if existing is None:
        xmat(filename, matrix, roles, good, nfull, starts)
    full = np.full((nfull, y.shape[1]), np.nan)
    full[good] = y
    volume(f'{tag}.nii', full)
    penalized = np.flatnonzero(roles)
    trials = np.column_stack([penalized, np.arange(len(penalized))+100,
                              conditions, trial_runs])
    # Shuffled input rows must not change column/beta identity.
    np.savetxt(f'{tag}.trials.1D', trials[::-1], fmt='%d',
               header='column trial condition run')
    args = ['-input', f'{tag}.nii', '-matrix', filename,
            '-fractions', ','.join(map(str, fractions)),
            '-trial_info', f'{tag}.trials.1D', '-save_cv', '-save_nuisance',
            '-save_fitted', '-save_residuals', '-write_manifest']
    run([BINARY, *args, '-prefix', tag, '-batch_size', 2, '-jobs', 2])
    loss, winners, pairs = repeat_oracle(matrix, y, roles, np.asarray(row_runs),
        np.asarray(trial_runs), np.asarray(conditions), fractions)
    dg = read(tag, 'diagnostics')
    assert_allclose(dg[10:], loss[:-1], rtol=4e-5, atol=4e-6)
    assert_allclose(dg[8], loss[-1], rtol=4e-5, atol=4e-6)
    assert_allclose(dg[9], winners, atol=0)
    assert_allclose(dg[0], np.asarray(fractions)[winners], atol=1e-7)
    assert_allclose(dg[7], loss[winners, np.arange(y.shape[1])],
                    rtol=4e-5, atol=4e-6)
    for v, win in enumerate(winners):
        fraction = fractions[win]
        expected = oracle(matrix, y[:, v:v+1], roles, fraction)
        for kind, values in zip(['betas', 'nuisance', 'fitted', 'residuals'], expected):
            assert_allclose(read(tag, kind)[:, v], values[:, 0], rtol=4e-5, atol=4e-6)
    info = nib.load(tag+'+orig.HEAD').header.info
    assert info['PICARD_TRIAL_INFO'] == ';'.join(' '.join(map(str,t)) for t in trials)
    assert 'fold-local' in info['PICARD_TUNING']
    assert 'Picard_CV_MSE' in info['BRICK_LABS']
    manifest = list(csv.DictReader(open(tag+'.betas.tsv'), delimiter='\t'))
    assert [int(t['Trial']) for t in manifest] == trials[:, 1].tolist()
    assert [int(t['Run']) for t in manifest] == trials[:, 3].tolist()
    run([BINARY, *args, '-prefix', tag+'_serial', '-batch_size', y.shape[1], '-jobs', 1])
    for kind in ['betas','nuisance','diagnostics','fitted','residuals']:
        assert np.array_equal(read(tag,kind),read(tag+'_serial',kind))
    return args, loss, pairs


def repeated_condition_tests(rng):
    ref = json.loads(Path(__file__).with_name('picard_glmsingle_reference.json').read_text())
    betas = np.asarray(ref['betas'])
    trial_runs = np.asarray(ref['run_ids'])
    conditions = np.asarray(ref['conditions'])
    p = len(betas)
    # Each trial occupies a +1/-1 pair, orthogonal to every other trial and
    # the three run intercepts. Distinct run means exercise nuisance recovery.
    x = np.zeros((2*p, 3+p))
    for j in range(p):
        x[2*j:2*j+2, trial_runs[j]] = 1
        x[2*j:2*j+2, 3+j] = [1,-1]
    roles = [0]*3+[1]*p
    y = x[:,3:]@betas + x[:,:3]@np.array([[10,20,30],[40,20,10],[-20,0,30]])
    args, losses, pairs = repeat_case('repeat_glm',x,y,roles,
        np.repeat(trial_runs,2),trial_runs,conditions,ref['fractions'])
    # Upstream z-scores with one session and returns SSE. This unit conversion
    # does not change the winner. No upstream code is needed at test time.
    converted = losses[:-1]*pairs/np.var(betas,axis=0,ddof=1)
    assert_allclose(converted,ref['badness'],rtol=3e-8,atol=3e-8)
    assert_allclose(read('repeat_glm','diagnostics')[9],ref['winners'])

    # Real 3dDeconvolve IM columns, uneven runs, censoring, overlapping HRFs,
    # an unequal-repeat condition and a response with no signal.
    Path('repeat_a.1D').write_text('1 5\n2\n1\n')
    Path('repeat_b.1D').write_text('3\n4\n3\n')
    n = 54; starts = [0,16,35]
    keep = np.ones(n,int); keep[[0,6,20,36]]=0
    good = np.flatnonzero(keep)
    np.savetxt('repeat_censor.1D',keep,fmt='%d')
    run([DECONVOLVE,'-nodata',n,1,'-polort',1,'-concat','1D: 0 16 35','-local_times',
         '-num_stimts',2,'-stim_times_IM',1,'repeat_a.1D','GAM',
         '-stim_label',1,'A','-stim_times_IM',2,'repeat_b.1D','GAM',
         '-stim_label',2,'B','-censor','repeat_censor.1D',
         '-x1D','repeat_real.xmat.1D','-x1D_stop'])
    x = np.loadtxt('repeat_real.xmat.1D').astype(np.float32).astype(float)
    assert x.shape == (len(good),13)
    roles = [0]*6+[1]*7
    row_runs = np.searchsorted(starts,good,side='right')-1
    trial_runs = np.array([0,0,1,2,0,1,2])
    conditions = np.array([0,0,0,0,1,1,1])
    y = x@rng.normal(size=(13,5))+rng.normal(size=(len(x),5))*.8
    y[:,4] = 0
    args, _, _ = repeat_case('repeat_real',x,y,roles,row_runs,trial_runs,
        conditions,[.2,.6,1,0],good,n,starts,existing='repeat_real.xmat.1D')
    assert read('repeat_real','diagnostics')[0,4] == 1
    # Perturb one held-out run strongly. The independent oracle computes all
    # training fits from the other rows; this catches full-data norm leakage.
    perturbed=y.copy(); perturbed[row_runs==1] *= -9
    repeat_case('repeat_perturbed',x,perturbed,roles,row_runs,trial_runs,
        conditions,[.15,.55,.9],good,n,starts)

    base = ['-input','repeat_real.nii','-matrix','repeat_real.xmat.1D']
    grid = ['-fractions','1,.5,0','-trial_info','repeat_real.trials.1D']
    failure_case([*base,*grid,'-fraction','.5'], 'either')
    failure_case([*base,'-fractions','1,,.5','-trial_info','repeat_real.trials.1D'], 'comma-separated')
    failure_case([*base,'-fractions','1,1','-trial_info','repeat_real.trials.1D'], 'duplicate')
    failure_case([*base,'-fractions','1,.5'], 'trial_info')
    failure_case([*base,*grid,'-memory_mb','1'], 'exceed')
    failure_case([*base,'-fraction','.5','-save_cv'], 'requires tuning')
    trial = np.loadtxt('repeat_real.trials.1D',dtype=int)
    for tag, column, values, message in [
        ('no_repeats',2,np.arange(len(trial)),'repeated'),
        ('same_trial',1,np.zeros(len(trial),int),'unique'),
        ('wrong_run',3,np.zeros(len(trial),int),'declared run'),
        ('wrong_column',0,np.zeros(len(trial),int),'trial_info')]:
        invalid=trial.copy(); invalid[:,column]=values
        np.savetxt(tag+'.1D',invalid,fmt='%d')
        failure_case([*base,'-fractions','1,.5','-trial_info',tag+'.1D'],message)
    for contents in ['6 1 0\n','6 1 0 0 garbage\n','2147483648 1 0 0\n']:
        Path('malformed.1D').write_text(contents)
        failure_case([*base,'-fractions','1,.5','-trial_info','malformed.1D'],'trial_info')
    # A nuisance term shared across runs would leak held-out information.
    shared=x.copy(); shared[:,0]=1
    xmat('repeat_shared.xmat.1D',shared,roles,good,n,starts)
    failure_case(['-input','repeat_real.nii','-matrix','repeat_shared.xmat.1D',*grid],
                 'run-local nuisance')
    validate=run([BINARY,'-matrix','repeat_real.xmat.1D',*grid,'-validate_only'])
    assert 'RepeatTuning: Folds=3 Candidates=3' in validate
    # Default tuned output is still just one bucket pair; score paths opt in.
    run([BINARY,*base,*grid,'-prefix','repeat_compact'])
    assert len(list(Path('.').glob('repeat_compact*')))==2
    assert read('repeat_compact','diagnostics').shape[0]==10
    volume('repeat_mask.nii',[[1,1,1,1,0]])
    run([BINARY,*base,*grid,'-prefix','repeat_masked','-mask','repeat_mask.nii'])
    assert np.all(read('repeat_masked','diagnostics')[7:,4]==-1)
    # Non-identifiable trials and redundant nuisance columns still use the
    # documented minimum-norm reference, independently in each training fold.
    deficient=x.copy(); deficient[:,7]=deficient[:,6]
    deficient=np.column_stack([deficient,deficient[:,0]])
    repeat_case('repeat_rank',deficient,y,roles+[0],row_runs,trial_runs,
        conditions,[1,.5,0],good,n,starts)

    # Grid order cannot change the selected fraction; index metadata changes.
    permuted=['-fractions','0,.5,1','-trial_info','repeat_real.trials.1D']
    run([BINARY,*base,*permuted,'-prefix','repeat_reordered'])
    assert_allclose(read('repeat_reordered','betas'),read('repeat_compact','betas'),atol=0)
    assert_allclose(read('repeat_reordered','diagnostics')[0],
                    read('repeat_compact','diagnostics')[0],atol=0)
    run([BINARY,*base,*grid,'-prefix','repeat_export.nii.gz'])
    exported=np.asarray(nib.load('repeat_export.nii.gz').dataobj).reshape(5,-1).T
    assert_allclose(exported[:7],read('repeat_compact','betas'),atol=0)
    Path('repeat_blocked').write_text('keep me')
    failure_case([*base,*grid,'-fitts','repeat_blocked/fitts'],prefix='repeat_rollback')
    assert Path('repeat_blocked').read_text()=='keep me'
    # Entirely censored/zero-support trials must not silently enter tuning.
    absent=x.copy(); absent[:,6]=0
    xmat('repeat_absent.xmat.1D',absent,roles,good,n,starts)
    failure_case(['-input','repeat_real.nii','-matrix','repeat_absent.xmat.1D',*grid],
                 'nonzero support')
    print('PASS repeated-condition fold oracle, pinned GLMsingle scoring, '
          'real IM/censor maps, final refits, failures, batches and threads')

def prediction_folds(good, starts, block=0, gap=0, groups=None):
    """Independent explicit index oracle in ORIGINAL acquisition coordinates."""
    runs=np.searchsorted(starts,good,side='right')-1
    if groups is not None:
        keys=list(groups)
    elif block:
        keys=[(int(r),(int(row)-starts[r])//block) for row,r in zip(good,runs)]
    else:
        keys=list(runs)
    unique=list(dict.fromkeys(keys))
    result=[]
    for key in unique:
        test=np.array([i for i,k in enumerate(keys) if k==key])
        train=np.array([i for i,k in enumerate(keys) if k!=key and not any(
            runs[i]==runs[j] and abs(int(good[i])-int(good[j]))<=gap for j in test)])
        assert len(train)>0
        result.append((train,test))
    return result


def prediction_select(loss, weights, fractions, one_se):
    """Loss shape [fold,candidate,response]; weighted-fold one-SE heuristic."""
    mean=np.einsum('f,fcv->cv',weights,loss)
    result=[]
    ordered=sorted(range(len(fractions)),key=lambda c:-fractions[c])
    for v in range(mean.shape[1]):
        best=min(ordered,key=lambda c:mean[c,v])
        se=np.sqrt(np.sum(weights*(loss[:,best,v]-mean[best,v])**2)*
                   np.sum(weights**2)/(1-np.sum(weights**2)))
        if one_se:
            eligible=np.flatnonzero(mean[:,v]<=mean[best,v]+se)
            best=min(eligible,key=lambda c:fractions[c])
        result.append(best)
    return np.array(result),mean


def prediction_oracle(x,y,roles,folds,fractions,project=False,fit=oracle):
    roles=np.asarray(roles); losses=[]; baselines=[]
    for train,test in folds:
        z=x[test][:,roles==0]; xt=x[test][:,roles==1]
        def residualize(values):
            if not project: return values.copy()
            r=values-z@np.linalg.lstsq(z,values,rcond=1e-10)[0]
            tiny=np.linalg.norm(r,axis=0)<=64*np.finfo(float).eps*max(z.shape)*np.linalg.norm(values,axis=0)
            r[:,tiny]=0
            return r
        yt=residualize(y[test]); xp=residualize(xt)
        baseline=residualize(np.ones((len(test),1))*y[train].mean(axis=0))
        baselines.append(np.mean((yt-baseline)**2,axis=0))
        fold=[]
        for fraction in fractions:
            beta,gamma,*_=fit(x[train],y[train],roles,fraction)
            fitted=xp@beta
            if not project: fitted+=z@gamma
            fold.append(np.mean((yt-fitted)**2,axis=0))
        losses.append(fold)
    return np.asarray(losses),np.asarray(baselines)


def prediction_case(tag,x,y,roles,good,starts,fractions,project=False,
                    one_se=False,block=0,gap=0,groups=None,global_mask=None,
                    existing=None):
    x=np.asarray(x,np.float32).astype(float); y=np.asarray(y,np.float32).astype(float)
    nfull=int(good[-1])+1
    filename=existing or tag+'.xmat.1D'
    if not existing: xmat(filename,x,roles,good,nfull,starts)
    full=np.full((nfull,y.shape[1]),np.nan); full[good]=y
    volume(tag+'.nii',full)
    args=['-input',tag+'.nii','-matrix',filename,'-cv_predict',
          '-fractions',','.join(map(str,fractions)),'-save_cv','-save_nuisance',
          '-save_fitted','-save_residuals','-cv_rule','1se' if one_se else 'min',
          '-cv_nuisance','project' if project else 'train','-cv_gap',gap]
    if block: args+=['-cv_block_size',block]
    if groups is not None:
        np.savetxt(tag+'.folds.1D',np.c_[good,groups][::-1],fmt='%d')
        args+=['-cv_folds',tag+'.folds.1D']
    if global_mask is not None:
        volume(tag+'.mask.nii',[global_mask])
        args+=['-cv_global','-mask',tag+'.mask.nii']
    folds=prediction_folds(good,starts,block,gap,groups)
    losses,baselines=prediction_oracle(x,y,roles,folds,fractions,project)
    weights=np.array([len(te) for _,te in folds])/len(good)
    baseline=weights@baselines
    selected,mean=prediction_select(losses,weights,fractions,one_se)
    keep=np.ones(y.shape[1],bool) if global_mask is None else np.asarray(global_mask,bool)
    defined=keep & (baseline>0)
    if global_mask is not None:
        if defined.any():
            normalized=(losses[:,:,defined]/baseline[defined]).mean(axis=2,keepdims=True)
            chosen,_=prediction_select(normalized,weights,fractions,one_se)
            selected[:]=chosen[0]
        else: selected[:]=int(np.argmax(fractions))
    run([BINARY,*args,'-prefix',tag,'-batch_size',2,'-jobs',2])
    dg=read(tag,'diagnostics')
    assert_allclose(dg[12:,keep],mean[:,keep],rtol=4e-5,atol=4e-6)
    assert_allclose(dg[7,keep],mean[selected,np.arange(y.shape[1])][keep],rtol=4e-5,atol=4e-6)
    assert_allclose(dg[8,keep],baseline[keep],rtol=4e-5,atol=4e-6)
    assert_allclose(dg[10,keep],selected[keep],atol=0)
    assert_allclose(dg[0,keep],np.asarray(fractions)[selected][keep],atol=1e-7)
    assert_allclose(dg[9,defined],1-dg[7,defined]/baseline[defined],rtol=4e-5,atol=4e-6)
    assert np.all(dg[11,keep & (baseline==0)]==1)
    assert np.all(dg[9,keep & (baseline==0)]==0)
    assert np.all(dg[11,~keep]==-1)
    for v in np.flatnonzero(keep):
        expected=oracle(x,y[:,v:v+1],roles,fractions[selected[v]])
        for kind,values in zip(['betas','nuisance','fitted','residuals'],expected):
            assert_allclose(read(tag,kind)[:,v],values[:,0],rtol=4e-5,atol=4e-6)
    header=nib.load(tag+'+orig.HEAD').header.info
    membership=header['PICARD_PRED_FOLDS'].split(';')
    assert len(membership)==len(folds)
    for record,(tr,te) in zip(membership,folds):
        assert list(map(int,record.split('train:')[1].split(' test:')[0].split(',')))==good[tr].tolist()
        assert list(map(int,record.split('test:')[1].split(',')))==good[te].tolist()
    assert header['PICARD_PRED_NUISANCE']==('project' if project else 'train')
    if global_mask is not None:
        assert f'defined_voxels={defined.sum()};' in header['PICARD_PRED_GLOBAL']
    run([BINARY,*args,'-prefix',tag+'_serial','-batch_size',y.shape[1],'-jobs',1])
    for kind in ['betas','nuisance','diagnostics','fitted','residuals']:
        assert np.array_equal(read(tag,kind),read(tag+'_serial',kind)),(tag,kind)
    return args,dg


def prediction_tests(rng):
    n=47; starts=[0,12,31]; good=np.delete(np.arange(n),[0,4,17,23,33])
    t=np.arange(n); features=np.c_[np.sin(t*.51),np.cos(t*.21)+.3*np.sin(t*.51)]
    x=np.c_[np.ones(len(good)),features[good]]
    roles=[0,1,1]; runs=np.searchsorted(starts,good,side='right')-1
    y=x@rng.normal(size=(3,5))+rng.normal(size=(len(x),5))*.8
    y[:,1]=features[good,0]*np.choose(runs,[5,-5,5])+rng.normal(size=len(x))*.1
    y[:,3]*=1000; y[:,4]=0
    args,dg=prediction_case('predict_raw',x,y,roles,good,starts,[1,.6,.2])
    assert (dg[9,:4]<0).any(), 'negative R2 must be retained'
    shifted=y.copy(); shifted[runs==1]*=-7
    prediction_case('predict_perturb',x,shifted,roles,good,starts,[.2,.6,1])
    prediction_case('predict_1se',x,y,roles,good,starts,[.2,.6,1,0],one_se=True)
    prediction_case('predict_global',x,y,roles,good,starts,[1,.6,.2,0],
                    one_se=True,global_mask=[1,1,1,1,0])
    prediction_case('predict_undefined',x,y,roles,good,starts,[.2,.6,1],
                    one_se=True,global_mask=[0,0,0,0,1])
    prediction_case('predict_blocks',x,y,roles,good,starts,[1,.5,0],block=6,gap=2)
    groups=np.array([(row-starts[r])//5%3+10 for row,r in zip(good,runs)])
    prediction_case('predict_manifest',x,y,roles,good,starts,[1,.5,0],groups=groups,gap=1)
    # Correlated rank-deficient predictors and redundant shared nuisance.
    redundant=np.c_[x,x[:,1]*2,x[:,0]*3]
    prediction_case('predict_rank',redundant,y,[0,1,1,1,0],good,starts,[1,.5,0])

    # Unchanged Deconvolve supplies ordinary shared stimulus columns and
    # run-local baselines. Their held-out intercepts require explicit project
    # scoring; raw prediction must reject them rather than silently fit Y_test.
    np.savetxt('predict_a.1D',features[:,0]);np.savetxt('predict_b.1D',features[:,1])
    keep=np.zeros(n,int); keep[good]=1;np.savetxt('predict_censor.1D',keep,fmt='%d')
    run([DECONVOLVE,'-nodata',n,1,'-polort',1,'-concat','1D: 0 12 31',
         '-num_stimts',2,'-stim_file',1,'predict_a.1D','-stim_file',2,'predict_b.1D',
         '-censor','predict_censor.1D','-x1D','predict_real.xmat.1D','-x1D_stop'])
    real=np.loadtxt('predict_real.xmat.1D').astype(np.float32).astype(float)
    real_y=real@rng.normal(size=(8,5))+rng.normal(size=(len(real),5))*.6;real_y[:,4]=0
    prediction_case('predict_real',real,real_y,[0]*6+[1]*2,good,starts,[1,.5,0],
                    project=True,existing='predict_real.xmat.1D')
    base=['-input','predict_raw.nii','-matrix','predict_raw.xmat.1D']
    mode=['-cv_predict','-fractions','1,.5,0']
    failure_case([*base,*mode,'-cv_global'],'explicit -mask')
    failure_case([*base,'-fraction','.5','-cv_rule','1se'],'require -cv_predict')
    failure_case([*base,'-fraction','.5','-cv_gap','0'],'require -cv_predict')
    failure_case([*base,*mode,'-cv_gap','-1'],'nonnegative')
    failure_case([*base,*mode,'-cv_rule','mystery'],'min or 1se')
    failure_case([*base,*mode,'-cv_nuisance','mystery'],'train or project')
    failure_case([*base,*mode,'-cv_folds','predict_manifest.folds.1D','-cv_block_size','6'],'blocks OR')
    failure_case([*base,*mode,'-cv_folds','predict_manifest.folds.1D','-cv_gap','1000'],'training')
    failure_case([*base,*mode,'-memory_mb','1'],'exceed')
    failure_case(['-input','predict_real.nii','-matrix','predict_real.xmat.1D',*mode],
                 'not estimable')
    failure_case(['-input','repeat_real.nii','-matrix','repeat_real.xmat.1D',
                  *mode,'-cv_nuisance','project'],'training support')
    for content in ['1 10\n1 11\n','0 10\n','1 10 extra\n','1 -1\n']:
        Path('bad_folds.1D').write_text(content)
        failure_case([*base,*mode,'-cv_folds','bad_folds.1D'],'cv_folds needs')
    validation=run([BINARY,'-matrix','predict_raw.xmat.1D',*mode,'-validate_only'])
    assert 'PredictionCV: Folds=3 Candidates=3' in validation
    run([BINARY,*base,*mode,'-prefix','predict_compact'])
    assert len(list(Path('.').glob('predict_compact*')))==2
    assert read('predict_compact','diagnostics').shape[0]==12
    print('PASS shared-predictor CV oracle, raw/conditional nuisance, one-SE/global '
          'selection, original-time gaps, manifests, failures and scheduling')

def elastic_oracle(matrix, y, roles, lam, rho, standardize):
    """Independent sklearn solve after explicit NumPy nuisance projection.

    sklearn alpha == Picard lambda; fit_intercept=False is essential because
    the nuisance span, including any desired intercept, is supplied explicitly.
    """
    from sklearn.linear_model import ElasticNet
    roles = np.asarray(roles)
    z, x = matrix[:, roles == 0], matrix[:, roles == 1]
    xp = x-z@np.linalg.lstsq(z,x,rcond=1e-10)[0]
    yp = y-z@np.linalg.lstsq(z,y,rcond=1e-10)[0]
    for projected, original in [(xp,x),(yp,y)]:
        tiny = np.linalg.norm(projected,axis=0) <= (
            64*np.finfo(float).eps*max(z.shape)*np.linalg.norm(original,axis=0))
        projected[:,tiny] = 0
    scale = np.sqrt(np.mean(xp*xp,axis=0)) if standardize else np.ones(x.shape[1])
    scale[scale == 0] = 1
    effective = xp/scale
    theta = []
    for v in range(y.shape[1]):
        if rho == 0:
            theta.append(np.linalg.solve(effective.T@effective +
                                         len(y)*lam*np.eye(x.shape[1]),
                                         effective.T@yp[:,v]))
        elif not np.any(yp[:,v]):
            theta.append(np.zeros(x.shape[1]))
        else:
            model = ElasticNet(alpha=lam,l1_ratio=rho,fit_intercept=False,
                               tol=1e-12,max_iter=100000,selection='cyclic')
            theta.append(model.fit(effective,yp[:,v]).coef_)
    theta = np.asarray(theta).T
    beta = theta/scale[:,None]
    gamma = np.linalg.lstsq(z,y-x@beta,rcond=1e-10)[0]
    fitted = x@beta+z@gamma
    residual = y-fitted
    gradient = -(effective.T@(yp-effective@theta))/len(y)+lam*(1-rho)*theta
    violation = np.where(theta != 0,abs(gradient+lam*rho*np.sign(theta)),
                         np.maximum(abs(gradient)-lam*rho,0)).max(axis=0)
    objective = np.mean(residual*residual,axis=0)/2+lam*(
        rho*abs(theta).sum(axis=0)+(1-rho)*(theta*theta).sum(axis=0)/2)
    maximum = abs(effective.T@yp).max(axis=0)/(len(y)*rho) if rho else -np.ones(y.shape[1])
    return beta,gamma,fitted,residual,objective,maximum,violation


def check_elastic(tag,x,y,roles,lambdas,rho,standardize,keep=None):
    keep = np.ones(y.shape[1],bool) if keep is None else keep
    dg = read(tag,'diagnostics')
    for v in np.flatnonzero(keep):
        result = elastic_oracle(x,y[:,v:v+1],roles,lambdas[v],rho,standardize)
        for kind,expected in zip(['betas','nuisance','fitted','residuals'],result):
            assert_allclose(read(tag,kind)[:,v],expected[:,0],rtol=3e-5,atol=3e-6)
        assert_allclose(dg[5,v],result[4][0],rtol=3e-5,atol=3e-6)
        assert_allclose(dg[6,v],result[5][0],rtol=3e-5,atol=3e-6)
        assert result[6][0] < 2e-8
    assert_allclose(dg[0,keep],np.asarray(lambdas)[keep],rtol=1e-6)
    assert_allclose(dg[1,keep],rho,rtol=1e-6)
    assert np.all(dg[4,keep] <= 1.01e-9)
    assert np.all(dg[9,keep] == 0)
    assert np.all(dg[9,~keep] == -1)
    assert np.array_equal(dg[2,keep],np.count_nonzero(read(tag,'betas')[:,keep],axis=0))
    hdr = nib.load(tag+'+orig.HEAD').header.info
    assert 'phase4-v1' in hdr['PICARD_MODEL']
    assert 'Picard_Fraction' not in hdr['BRICK_LABS']
    assert 'Picard_L1Ratio' in hdr['BRICK_LABS']


def elastic_tests(rng):
    # Unequal runs, removed rows, different feature scales across runs, a
    # redundant nuisance column and a zero response exercise the full adapter.
    starts = np.array([0,20,43]); good = np.delete(np.arange(66),[3,24,48])
    run_ids = np.searchsorted(starts,good,side='right')-1
    x = rng.normal(size=(len(good),5))*[.3,1,2,4,1]
    x[:,1] += .6*x[:,2]
    x[run_ids == 1,0] *= 9
    z = np.c_[np.ones(len(good)),good/66,np.ones(len(good))*2]
    matrix = xmat('elastic.xmat.1D',np.c_[x,z],[1]*5+[0]*3,good,66,starts)
    roles = np.array([1]*5+[0]*3)
    y = (matrix[:,:5]@rng.normal(size=(5,4)) +
         matrix[:,5:]@rng.normal(size=(3,4))+rng.normal(size=(len(good),4))*.5)
    y[:,3] = 0
    y = y.astype(np.float32).astype(float)
    full = np.full((66,4),np.nan); full[good] = y
    volume('elastic.nii',full)
    base = ['-input','elastic.nii','-matrix','elastic.xmat.1D',
            '-solver_tol','1e-9','-save_nuisance','-save_fitted','-save_residuals']
    for rho in [0,.4,1]:
        for scaled in [False,True]:
            tag = f'elastic_fixed_{rho}_{scaled}'
            flags = ['-penalty','lasso' if rho == 1 else 'enet','-l1_ratio',rho]
            if scaled: flags += ['-standardize']
            run([BINARY,*base,*flags,'-lambda',.12,'-prefix',tag])
            check_elastic(tag,matrix,y,roles,[.12]*4,rho,scaled)
    # Raw and conditional prediction, global/per-voxel selection, both rules.
    grid = np.array([1.,.2,.03])
    weights = np.array([sum(run_ids==r) for r in range(3)])/len(good)
    folds = prediction_folds(good,starts)
    for project,one_se,global_choice,rho in [
            (False,False,False,1), (False,True,True,.4),
            (True,False,True,1), (True,True,False,.4)]:
        tag = f'elastic_cv_{project}_{one_se}_{global_choice}'
        def fit(x,y,roles,lam):
            return elastic_oracle(x,y,roles,lam,rho,True)
        losses,baselines = prediction_oracle(matrix,y,roles,folds,grid,project,fit)
        # Existing fraction selector is parameterized by NEGATIVE lambda so
        # its weak/strong ordering agrees; it does not perform any fitting.
        selected,mean = prediction_select(losses,weights,-grid,one_se)
        baseline = weights@baselines
        args = [*base,'-penalty','enet','-l1_ratio',rho,'-standardize',
                '-cv_predict','-lambdas',','.join(map(str,grid)),'-save_cv',
                '-cv_rule','1se' if one_se else 'min',
                '-cv_nuisance','project' if project else 'train']
        if global_choice:
            volume(tag+'.mask.nii',[[1,1,1,1]])
            args += ['-cv_global','-mask',tag+'.mask.nii']
            valid = baseline > 0
            normalized = (losses[:,:,valid]/baseline[valid]).mean(axis=2,keepdims=True)
            selected[:] = prediction_select(normalized,weights,-grid,one_se)[0][0]
        run([BINARY,*args,'-prefix',tag,'-jobs',2,'-batch_size',2])
        dg = read(tag,'diagnostics')
        assert_allclose(dg[15:],mean,rtol=4e-5,atol=4e-6)
        assert_allclose(dg[10],mean[selected,np.arange(4)],rtol=4e-5,atol=4e-6)
        assert_allclose(dg[11],baseline,rtol=4e-5,atol=4e-6)
        assert np.array_equal(dg[13],selected)
        check_elastic(tag,matrix,y,roles,grid[selected],rho,True)
        hdr = nib.load(tag+'+orig.HEAD').header.info
        assert 'PICARD_LAMBDAS' in hdr and 'PICARD_FRACTIONS' not in hdr
        run([BINARY,*args,'-prefix',tag+'_serial','-jobs',1,'-batch_size',4])
        for kind in ['betas','nuisance','fitted','residuals','diagnostics']:
            assert np.array_equal(read(tag,kind),read(tag+'_serial',kind))
    # Trial-unique columns use beta replication, never time-series prediction.
    trials = np.zeros((36,6)); nuisance = np.zeros((36,3))
    for r in range(3):
        trials[r*12:(r+1)*12,r*2:r*2+2] = rng.normal(size=(12,2))
        nuisance[r*12:(r+1)*12,r] = 1
    trial_matrix = xmat('elastic_trials.xmat.1D',np.c_[trials,nuisance],
                        [1]*6+[0]*3,runs=[0,12,24])
    trial_y = (trial_matrix[:,:6]@np.array([[1,2],[.5,-1]]*3)+
               rng.normal(size=(36,2))*.8).astype(np.float32).astype(float)
    volume('elastic_trials.nii',trial_y)
    np.savetxt('elastic_trials.1D',np.c_[np.arange(6),np.arange(6),
               np.tile([0,1],3),np.repeat(np.arange(3),2)],fmt='%d')
    target = oracle(trial_matrix,trial_y,[1]*6+[0]*3,1)[0]
    loss = np.zeros((len(grid),2))
    for r in range(3):
        train_rows = np.repeat(np.arange(3),12) != r
        train_cols = np.r_[np.repeat(np.arange(3),2) != r,np.arange(3) != r]
        for c,lam in enumerate(grid):
            beta = elastic_oracle(trial_matrix[train_rows][:,train_cols],
                trial_y[train_rows],[1]*4+[0]*2,lam,.6,True)[0]
            for j,original in enumerate(np.flatnonzero(train_cols[:6])):
                loss[c] += (beta[j]-target[2*r+original%2])**2/12
    winners = np.array([min(range(len(grid)),key=lambda c:(loss[c,v],grid[c]))
                        for v in range(2)])
    run([BINARY,'-input','elastic_trials.nii','-matrix','elastic_trials.xmat.1D',
         '-penalty','enet','-l1_ratio',.6,'-standardize','-solver_tol','1e-9',
         '-lambdas',','.join(map(str,grid)),'-trial_info','elastic_trials.1D',
         '-save_cv','-save_nuisance','-save_fitted','-save_residuals',
         '-prefix','elastic_repeat'])
    dg = read('elastic_repeat','diagnostics')
    assert_allclose(dg[13:],loss,rtol=4e-5,atol=4e-6)
    assert np.array_equal(dg[12],winners)
    check_elastic('elastic_repeat',trial_matrix,trial_y,[1]*6+[0]*3,
                  grid[winners],.6,True)
    for options in [
        ['-penalty','lasso','-lambda',0],
        ['-penalty','enet','-lambda',.1,'-l1_ratio',2],
        ['-penalty','lasso','-lambda',.1,'-l1_ratio',.5],
        ['-penalty','enet','-fraction',.5],
        ['-penalty','ridge','-lambda',.1],
        ['-penalty','enet','-lambda',.1,'-max_iterations',1,'-solver_tol','1e-30'],
        ['-penalty','enet','-lambdas','.1,1','-cv_predict'],
        ['-penalty','enet','-lambdas','1,1','-cv_predict'],
        ['-penalty','enet','-lambdas','1,0','-cv_predict']]:
        failure_case([*base,*options])
    # Wide, rank-deficient lasso: compare predictions/objectives because
    # duplicate columns admit different equally optimal coefficient vectors.
    wide = rng.normal(size=(8,10)); wide[:,1] = wide[:,0]; wide[:,2] = 0
    wide = xmat('elastic_wide.xmat.1D',np.c_[wide,np.ones(8)],[1]*10+[0])
    wy = rng.normal(size=(8,2)).astype(np.float32).astype(float)
    volume('elastic_wide.nii',wy)
    run([BINARY,'-input','elastic_wide.nii','-matrix','elastic_wide.xmat.1D',
         '-penalty','lasso','-lambda',.1,'-solver_tol','1e-9','-standardize',
         '-save_fitted','-prefix','elastic_wide'])
    wr = elastic_oracle(wide,wy,[1]*10+[0],.1,1,True)
    assert_allclose(read('elastic_wide','fitted'),wr[2],atol=5e-6,rtol=5e-5)
    assert_allclose(read('elastic_wide','diagnostics')[5],wr[4],atol=5e-6)
    # All undefined baselines fall back to the weakest lambda. Explicit
    # temporal folds and gaps still use original time, including censored rows.
    volume('elastic_zero.nii',np.zeros((66,2)))
    volume('elastic_zero.mask.nii',[[1,1]])
    run([BINARY,'-input','elastic_zero.nii','-matrix','elastic.xmat.1D',
         '-penalty','enet','-lambdas','2,.4,.01','-cv_predict','-cv_global',
         '-mask','elastic_zero.mask.nii','-cv_rule','1se','-cv_block_size',10,
         '-cv_gap',1,'-prefix','elastic_zero'])
    zd = read('elastic_zero','diagnostics')
    assert_allclose(zd[0],.01,rtol=1e-6)
    assert np.all(zd[13] == 2) and np.all(zd[14] == 1)
    # A real unchanged 3dDeconvolve design also works in the new fixed mode.
    # Its earlier fixture includes a censored row and run-local nuisance terms.
    real = Path('real.xmat.1D')
    assert real.exists()
    run([BINARY,'-input','input.nii','-mask','mask.nii','-matrix',real,'-penalty','lasso',
         '-lambda',.1,'-prefix','elastic_real'])
    print('PASS sklearn elastic-net/lasso, ridge endpoint, scaling, both CV modes, '
          'global/one-SE, batch/thread invariance, diagnostics and failure rollback')


def candidate_manifest(tag, matrices, roles, good, starts, ids=(11,4,9),
                       trial_runs=None, conditions=None, items=None):
    """Versioned, relative-path fixture with explicit per-candidate identities."""
    folder = Path(tag+' designs'); folder.mkdir()
    ncol = len(roles); labels = [f'identity{j}' for j in range(ncol)]
    entries = []
    arrays = []
    for c,(matrix,identifier) in enumerate(zip(matrices,ids)):
        name = f'design {identifier}.xmat.1D'
        arrays.append(xmat(folder/name,matrix,roles,good,int(good[-1])+1,starts,labels))
        table = folder/f'columns {identifier}.tsv'
        rows = ['column\tlabel\trole\ttrial\tcondition\trun\tbasis\titem']
        trial = 0
        for j,role in enumerate(roles):
            identity = ['-','-','-']
            if role and trial_runs is not None:
                identity = [str(trial),str(conditions[trial]),str(trial_runs[trial])]
                trial += 1
            rows.append('\t'.join([str(j),labels[j],
                'penalized' if role else 'nuisance',*identity,'0',
                items[j] if items else '-']))
        # Identity file order is deliberately unrelated to column order.
        table.write_text(rows[0]+'\n'+'\n'.join(rows[:0:-1])+'\n')
        entries.append('\t'.join([str(identifier),name,table.name,
                       'synthetic-reference','v1','0.1','already-convolved',
                       'shared-original-units']))
    manifest = folder/'candidates.tsv'
    manifest.write_text('# 3dPicard candidates 1\n'
        'id\tmatrix\tcolumns\tsource\tversion\tsampling\tduration\tnormalization\n'+
        '\n'.join(entries)+'\n')
    return manifest, arrays


def candidate_maps(tag):
    image = nib.load(tag+'+orig.HEAD'); header=image.header.info
    data=np.asarray(image.dataobj).reshape(image.shape[0],-1).T
    start=header['PICARD_DESIGN_DIAGNOSTIC_START']
    return data[start:],header


def candidate_tests(rng):
    # Two distinct shared-predictor designs, and an exact copy under a smaller
    # stable ID. Voxels prefer different designs; zero response is an exact tie.
    good=np.delete(np.arange(54),[3,26,40]); starts=[0,18,36]
    z=np.c_[np.ones(len(good)),np.linspace(-1,1,len(good))]
    first=rng.normal(size=(len(good),3)); second=rng.normal(size=first.shape)
    matrices=[np.c_[first,z],np.c_[second,z],np.c_[first,z]]
    roles=[1,1,1,0,0]; ids=np.array([11,4,9])
    manifest,matrices=candidate_manifest('candidate',matrices,roles,good,starts)
    y=np.c_[matrices[0][:,0]*2,matrices[1][:,1]*3,matrices[0][:,2]*1.5,
            np.zeros(len(good))]+z@np.array([[.5,1,-.4,0],[.1,.2,.1,0]])
    y[:,:3]+=rng.normal(size=(len(good),3))*.03
    y=y.astype(np.float32).astype(float)
    full=np.full((54,4),np.nan);full[good]=y;volume('candidate.nii',full)
    base=['-input','candidate.nii','-candidate_manifest',manifest,
          '-save_candidates','-save_nuisance','-save_fitted','-save_residuals',
          '-write_manifest']
    ols_loss=np.array([np.mean(oracle(x,y,roles,1)[3]**2,axis=0) for x in matrices])
    selected=np.array([min(range(3),key=lambda c:(ols_loss[c,v],ids[c]))
                       for v in range(4)])
    assert ids[selected].tolist()==[9,4,9,4]
    for fraction in [0,.6,1]:
        tag=f'candidate_hrf_{fraction}'
        args=[*base,'-candidate_policy','hrf','-fraction',fraction]
        run([BINARY,*args,'-prefix',tag,'-jobs',2,'-batch_size',3])
        maps,header=candidate_maps(tag)
        assert np.array_equal(maps[0],ids[selected])
        assert np.array_equal(maps[1],selected)
        assert_allclose(maps[2],ols_loss[selected,np.arange(4)],atol=2e-7)
        assert_allclose(maps[3:],ols_loss,rtol=4e-5,atol=4e-6)
        for v,c in enumerate(selected):
            expected=oracle(matrices[c],y[:,v:v+1],roles,fraction)
            for kind,value in zip(['betas','nuisance','fitted','residuals'],expected):
                assert_allclose(read(tag,kind)[:,v],value[:,0],rtol=4e-5,atol=4e-6)
        assert 'full-data' in header['PICARD_CANDIDATE_SELECTION']
        assert header['PICARD_MATRIX']=='MULTIPLE: see PICARD_CANDIDATES'
        assert str(manifest.resolve())==header['PICARD_CANDIDATE_MANIFEST']
        table=list(csv.DictReader(open(tag+'.betas.tsv'),delimiter='\t'))
        assert all('Xmat' not in row and 'XColumn' not in row for row in table)
        assert [int(row['Column']) for row in table]==[0,1,2]
        assert int(table[0]['DesignIDVolume'])==header['PICARD_DESIGN_DIAGNOSTIC_START']
        run([BINARY,*args,'-prefix',tag+'_serial','-jobs',1,'-batch_size',4])
        for kind in ['betas','nuisance','fitted','residuals','diagnostics']:
            assert np.array_equal(read(tag,kind),read(tag+'_serial',kind))
    # HRF choice is independent of the final penalty family and its strength.
    tag='candidate_hrf_lasso'
    run([BINARY,*base,'-candidate_policy','hrf','-penalty','lasso','-lambda',.1,
         '-standardize','-solver_tol','1e-9','-prefix',tag])
    assert np.array_equal(candidate_maps(tag)[0][0],ids[selected])
    for v,c in enumerate(selected):
        expected=elastic_oracle(matrices[c],y[:,v:v+1],roles,.1,1,True)
        assert_allclose(read(tag,'betas')[:,v],expected[0][:,0],rtol=4e-5,atol=4e-6)
    # Candidate ordering changes indices and score-map order, never stable IDs.
    lines=manifest.read_text().splitlines()
    reversed_manifest=manifest.with_name('reversed.tsv')
    reversed_manifest.write_text('\n'.join(lines[:2]+lines[:1:-1])+'\n')
    run([BINARY,'-input','candidate.nii','-candidate_manifest',reversed_manifest,
         '-candidate_policy','hrf','-fraction',.6,'-prefix','candidate_reordered'])
    assert np.array_equal(candidate_maps('candidate_reordered')[0][0],ids[selected])
    # Independently tune each design, then compare selected prediction losses.
    folds=prediction_folds(good,starts); weights=np.array([len(te) for _,te in folds])/len(good)
    for elastic,one_se,project in [(False,False,False),(False,True,True),
                                  (True,False,True),(True,True,False)]:
        grid=np.array([1,.6,.2]) if not elastic else np.array([.8,.2,.03])
        tag=f'candidate_predict_{elastic}_{one_se}_{project}'
        losses=[]; choices=[]; pooled=[]; baselines=[]
        def fit(x,y,roles,parameter):
            return (elastic_oracle(x,y,roles,parameter,.5,True) if elastic
                    else oracle(x,y,roles,parameter))
        for x in matrices:
            loss,baseline=prediction_oracle(x,y,roles,folds,grid,project,fit)
            choice,mean=prediction_select(loss,weights,-grid if elastic else grid,one_se)
            losses.append(mean[choice,np.arange(4)])
            choices.append(choice);pooled.append(mean);baselines.append(weights@baseline)
        losses=np.asarray(losses)
        winner=np.array([min(range(3),key=lambda c:(losses[c,v],ids[c]))
                         for v in range(4)])
        args=[*base,'-candidate_policy','predict','-cv_predict','-save_cv',
              '-cv_rule','1se' if one_se else 'min','-cv_nuisance',
              'project' if project else 'train']
        if elastic:args+=['-penalty','enet','-l1_ratio',.5,'-standardize',
                          '-solver_tol','1e-9','-lambdas',','.join(map(str,grid))]
        else:args+=['-fractions',','.join(map(str,grid))]
        run([BINARY,*args,'-prefix',tag,'-jobs',2,'-batch_size',2])
        maps,header=candidate_maps(tag); dg=read(tag,'diagnostics'); offset=10 if elastic else 7
        assert np.array_equal(maps[0],ids[winner])
        assert_allclose(maps[3:],losses,rtol=5e-5,atol=5e-6)
        assert 'outer-test' in header['PICARD_CANDIDATE_SELECTION']
        for v,c in enumerate(winner):
            parameter=grid[choices[c][v]]
            result=fit(matrices[c],y[:,v:v+1],roles,parameter)
            for kind,value in zip(['betas','nuisance','fitted','residuals'],result):
                assert_allclose(read(tag,kind)[:,v],value[:,0],rtol=5e-5,atol=5e-6)
            assert_allclose(dg[offset+5:offset+5+len(grid),v],pooled[c][:,v],rtol=5e-5,atol=5e-6)
            assert dg[offset+3,v]==choices[c][v]
        run([BINARY,*args,'-prefix',tag+'_serial','-jobs',1,'-batch_size',4])
        assert np.array_equal(read(tag,'diagnostics'),read(tag+'_serial','diagnostics'))
    # Fixed HRF candidate identities also survive masking and optional NIfTI.
    volume('candidate.mask.nii',[[1,0,1,0]])
    run([BINARY,*base,'-candidate_policy','hrf','-fraction',.6,
         '-mask','candidate.mask.nii','-prefix','candidate_mask'])
    assert np.all(candidate_maps('candidate_mask')[0][:,[1,3]]==-1)
    run([BINARY,'-input','candidate.nii','-candidate_manifest',manifest,
         '-candidate_policy','hrf','-fraction',.6,'-prefix','candidate.nii.gz'])
    ni=np.asarray(nib.load('candidate.nii.gz').dataobj).reshape(4,-1).T
    assert_allclose(ni[:3],read('candidate_hrf_0.6','betas'),atol=1e-6)
    # Fail before output on malformed identities, row/censor/run/TR mismatches,
    # nuisance changes, rank differences, mixed policy flags and low budgets.
    validate=['-candidate_manifest',manifest,'-candidate_policy','hrf','-validate_only']
    assert 'CandidateDesigns=3' in run([BINARY,*validate])
    failbase=['-input','candidate.nii','-candidate_manifest',manifest,
              '-candidate_policy','hrf','-fraction',.6]
    failure_case(failbase+['-matrix','candidate designs/design 11.xmat.1D'])
    failure_case(failbase+['-cv_global','-mask','candidate.mask.nii'])
    failure_case(failbase+['-memory_mb','.001'])
    failure_case(['-matrix','candidate designs/design 11.xmat.1D','-validate_only','-save_candidates'])
    # Mutations restore original bytes between independent rejection cases.
    target=manifest.parent/'columns 4.tsv'; original=target.read_text()
    for changed in [original.replace('identity0','different'),
                    original.replace('penalized','nuisance',1),
                    '\n'.join(original.splitlines()[:-1])+'\n',
                    original+original.splitlines()[1]+'\n',
                    original.replace('\t0\t-','\t1\t-',1)]:
        target.write_text(changed);run([BINARY,*validate],ok=False)
    target.write_text(original)
    matrix_file=manifest.parent/'design 4.xmat.1D'; original=matrix_file.read_text()
    for changed in [original.replace('0,18,36','0,19,36'),
                    original.replace('ColumnLabels = "','RowTR = "2"\n# ColumnLabels = "'),
                    original.replace('0,1,2,4','0,1,3,4')]:
        assert changed!=original
        matrix_file.write_text(changed);run([BINARY,*validate],ok=False)
    matrix_file.write_text(original)
    for kind in ['nuisance','rank']:
        bad=matrices[1].copy()
        if kind=='nuisance':bad[:,3]*=2
        else:bad[:,0]=0
        xmat(matrix_file,bad,roles,good,54,starts,[f'identity{j}' for j in range(5)])
        run([BINARY,*validate],ok=False)
    matrix_file.write_text(original)
    original=manifest.read_text()
    for changed in [original.replace('candidates 1','candidates 2'),
                    original.replace('\n4\t','\n11\t'),
                    original.replace('already-convolved','different',1),
                    original.replace('\n4\t','\n16777216\t')]:
        manifest.write_text(changed);run([BINARY,*validate],ok=False)
    manifest.write_text(original)
    # Single-trial HRF candidates retain one event identity per column. Choose
    # the HRF by OLS first, then tune repeated-condition betas on that design.
    from scipy.stats import gamma as gamma_distribution
    run_starts=[0,32,64]; timeline=np.arange(96); trial_runs=np.repeat(np.arange(3),2)
    conditions=np.tile([0,1],3); trial_roles=[1]*6+[0]*3
    generated=[]
    for shape in [4,8,4]:
        hrf=gamma_distribution.pdf(np.arange(24),shape);hrf/=hrf.max()
        design=np.zeros((96,9))
        for run_id in range(3):
            design[run_id*32:(run_id+1)*32,6+run_id]=1
            for trial,onset in enumerate([2,9]):
                event=np.zeros(32);event[onset]=1
                design[run_id*32:(run_id+1)*32,run_id*2+trial]=np.convolve(event,hrf)[:32]
        generated.append(design)
    tmanifest,generated=candidate_manifest('trial_candidate',generated,trial_roles,
        timeline,run_starts,trial_runs=trial_runs,conditions=conditions)
    betas=np.array([1,2,1,2,1,2.])
    ty=np.c_[generated[0][:,:6]@betas,generated[1][:,:6]@betas,
             generated[0][:,:6]@(-betas),np.zeros(96)]
    ty[:,:3]+=rng.normal(size=(96,3))*.02
    ty=ty.astype(np.float32).astype(float);volume('candidate_trials.nii',ty)
    np.savetxt('candidate_trials.1D',np.c_[np.arange(6),np.arange(6),conditions,trial_runs],fmt='%d')
    hrf_loss=np.array([np.mean(oracle(x,ty,trial_roles,1)[3]**2,axis=0) for x in generated])
    hrf_choice=np.array([min(range(3),key=lambda c:(hrf_loss[c,v],ids[c])) for v in range(4)])
    assert ids[hrf_choice].tolist()==[9,4,9,4]
    fractions=np.array([1,.7,.3])
    repeat_results=[repeat_oracle(x,ty,trial_roles,np.repeat(np.arange(3),32),
                    trial_runs,conditions,fractions) for x in generated]
    targs=['-input','candidate_trials.nii','-candidate_manifest',tmanifest,
           '-candidate_policy','hrf','-fractions','1,.7,.3','-trial_info','candidate_trials.1D',
           '-save_cv','-save_candidates','-save_nuisance','-save_fitted','-save_residuals',
           '-write_manifest']
    run([BINARY,*targs,'-prefix','candidate_repeat','-jobs',2,'-batch_size',3])
    maps,header=candidate_maps('candidate_repeat');dg=read('candidate_repeat','diagnostics')
    assert np.array_equal(maps[0],ids[hrf_choice])
    assert 'CONDITIONAL' in header['PICARD_CANDIDATE_SELECTION']
    for v,c in enumerate(hrf_choice):
        loss,winners,_=repeat_results[c]
        assert_allclose(dg[10:13,v],loss[:3,v],rtol=5e-5,atol=5e-6)
        assert dg[9,v]==winners[v]
        result=oracle(generated[c],ty[:,v:v+1],trial_roles,fractions[winners[v]])
        for kind,value in zip(['betas','nuisance','fitted','residuals'],result):
            assert_allclose(read('candidate_repeat',kind)[:,v],value[:,0],rtol=5e-5,atol=5e-6)
    table=list(csv.DictReader(open('candidate_repeat.betas.tsv'),delimiter='\t'))
    assert [int(row['Trial']) for row in table]==list(range(6))
    assert [int(row['Condition']) for row in table]==conditions.tolist()
    run([BINARY,*targs,'-prefix','candidate_repeat_serial','-jobs',1,'-batch_size',4])
    assert np.array_equal(read('candidate_repeat','diagnostics'),
                          read('candidate_repeat_serial','diagnostics'))
    # Alter one trial identity while preserving dimensions/labels: rejected.
    table_file=tmanifest.parent/'columns 4.tsv'; original=table_file.read_text()
    table_file.write_text(original.replace('penalized\t0\t0\t0',
                                           'penalized\t10\t0\t0'))
    run([BINARY,*targs,'-validate_only'],ok=False)
    table_file.write_text(original)
    print('PASS candidate HRF/prediction independent oracles, stable IDs, '
          'column identities, provenance, masks, ties, policies and failures')


def phase6_whitener(rows, runs, rho):
    """Independent dense matrix; never call production whitening code."""
    n=len(rows); w=np.eye(n)*np.sqrt(1-rho*rho)
    for i in range(1,n):
        if runs[i]==runs[i-1] and rows[i]==rows[i-1]+1:
            w[i,i]=1;w[i,i-1]=-rho
    return w


def phase6_fit(matrix,y,roles,rows,runs,rho,kind,parameter,groups=None,
               support=False,standardize=False):
    roles=np.asarray(roles);w=phase6_whitener(rows,runs,rho)
    wx,wy=w@matrix,w@y
    support_size=support_rank=None
    if kind=='banded':
        # Joint augmented LS is independent of Picard's projection and
        # primal/dual normal-equation implementation.
        diagonal=np.zeros(matrix.shape[1])
        diagonal[roles==1]=len(y)*np.asarray(parameter)[groups]
        coefficients=np.linalg.lstsq(np.r_[wx,np.diag(np.sqrt(diagonal))],
            np.r_[wy,np.zeros((matrix.shape[1],y.shape[1]))],rcond=1e-10)[0]
        beta,gamma=coefficients[roles==1],coefficients[roles==0]
    else:
        fitted=(oracle(wx,wy,roles,parameter) if kind=='ridge' else
                elastic_oracle(wx,wy,roles,parameter,.7,standardize))
        beta,gamma=fitted[:2]
        if support:
            selection=beta!=0
            support_size=selection.sum(axis=0);support_rank=[]
            beta=np.zeros_like(beta);gamma=np.zeros_like(gamma)
            for v in range(y.shape[1]):
                columns=roles==0
                columns[roles==1]=selection[:,v]
                coefficients=np.linalg.lstsq(wx[:,columns],wy[:,v],rcond=1e-10)[0]
                beta[selection[:,v],v]=coefficients[roles[columns]==1]
                gamma[:,v]=coefficients[roles[columns]==0]
                z=wx[:,roles==0];xp=wx[:,roles==1][:,selection[:,v]]
                xp=xp-z@np.linalg.lstsq(z,xp,rcond=1e-10)[0]
                support_rank.append(np.linalg.matrix_rank(xp) if xp.size else 0)
    fitted=matrix[:,roles==1]@beta+matrix[:,roles==0]@gamma
    return beta,gamma,fitted,y-fitted,support_size,support_rank


def phase6_check(prefix,expected):
    for kind,value in zip(['betas','nuisance','fitted','residuals'],expected):
        assert_allclose(read(prefix,kind),value,rtol=7e-5,atol=6e-6)


def phase6_tests(rng):
    good=np.delete(np.arange(63),[3,4,22,45]);starts=[0,21,42]
    runs=np.searchsorted(starts,good,side='right')-1
    roles=np.array([1,0,1,1,0,1]);groups=np.array([0,1,0,1])
    raw=rng.normal(size=(len(good),6));raw[:,1]=1;raw[:,4]=good/63
    raw[:,3]*=3;raw[:,5]*=.2
    x=xmat('p6.xmat.1D',raw,roles,good,63,starts)
    y=x@rng.normal(size=(6,5))+rng.normal(size=(len(good),5))*.4;y[:,4]=0
    y=y.astype(np.float32).astype(float)
    full=np.zeros((63,5));full[good]=y;volume('p6.nii',full)
    np.savetxt('p6.groups',[[0,-1,1,0,-1,1]],fmt='%d')
    vectors=np.array([[.03,.7],[.6,.02],[1.4,2]])
    np.savetxt('p6.penalties',vectors)
    np.savetxt('p6.fixed',vectors[:1])
    base=['-input','p6.nii','-matrix','p6.xmat.1D','-save_nuisance',
          '-save_fitted','-save_residuals']
    band=['-penalty','banded','-band_groups','p6.groups']
    for rho in [0,.55,-.4]:
        tag=f'p6_fixed_{rho}'
        flags=[*band,'-band_penalties','p6.fixed','-prewhite_ar1',rho]
        run([BINARY,*base,*flags,'-prefix',tag,'-batch_size',2,'-jobs',2])
        expected=phase6_fit(x,y,roles,good,runs,rho,'banded',vectors[0],groups)
        phase6_check(tag,expected)
        dg=read(tag,'diagnostics');w=phase6_whitener(good,runs,rho)
        objective=np.mean((w@expected[3])**2,axis=0)/2+np.sum(
            vectors[0][groups,None]*expected[0]**2,axis=0)/2
        assert_allclose(dg[1],objective,rtol=5e-5,atol=3e-6)
        assert np.max(dg[2])<1e-8 and np.all(dg[6]==0)
        assert_allclose(dg[4:6],np.repeat(vectors[0][:,None],5,axis=1))
        hdr=nib.load(tag+'+orig.HEAD').header.info
        assert 'PICARD_BAND_GROUPS' in hdr and 'PICARD_BAND_PENALTIES' in hdr
        assert 'FixedAR1=' in hdr['PICARD_MODEL']
        run([BINARY,*base,*flags,'-prefix',tag+'_serial','-batch_size',5,'-jobs',1])
        for kind in ['betas','nuisance','diagnostics','fitted','residuals']:
            assert np.array_equal(read(tag,kind),read(tag+'_serial',kind))
    # Explicit dual-system fixture: wide, duplicate predictors, no nuisance.
    wide=rng.normal(size=(5,9));wide[:,8]=wide[:,0]
    wide=xmat('p6wide.xmat.1D',wide,[1]*9)
    wy=rng.normal(size=(5,3)).astype(np.float32).astype(float)
    volume('p6wide.nii',wy);wg=np.arange(9)%2
    np.savetxt('p6wide.groups',wg,fmt='%d')
    run([BINARY,'-matrix','p6wide.xmat.1D','-input','p6wide.nii','-penalty','banded',
         '-band_groups','p6wide.groups','-band_penalties','p6.fixed',
         '-prefix','p6wide','-save_fitted','-save_residuals'])
    phase6_check('p6wide',phase6_fit(wide,wy,[1]*9,np.arange(5),np.zeros(5),
                                  0,'banded',vectors[0],wg))
    # Fixed fractional whitening and sparse support refits use independent
    # dense whitening, sklearn selection, and joint minimum-norm LS recovery.
    for kind,parameter,support in [('ridge',.6,False),('enet',.08,False),
                                    ('enet',.08,True),('enet',1e5,True)]:
        tag=f'p6_{kind}_{parameter}_{support}'
        flags=['-fraction',parameter] if kind=='ridge' else [
            '-penalty','enet','-lambda',parameter,'-l1_ratio',.7,
            '-standardize','-solver_tol',1e-10]
        if support:flags+=['-support_refit']
        run([BINARY,*base,*flags,'-prewhite_ar1',.55,'-prefix',tag])
        expected=phase6_fit(x,y,roles,good,runs,.55,kind,parameter,
                           support=support,standardize=True)
        phase6_check(tag,expected)
        if support:
            dg=read(tag,'diagnostics')
            assert np.array_equal(dg[9],expected[4])
            assert np.array_equal(dg[10],expected[5]) and np.all(dg[11]==0)
    # Duplicate selected columns test minimum-norm support recovery instead
    # of assuming that a selected support is automatically full rank.
    t=np.arange(6,dtype=float)-2
    duplicate=xmat('p6support_rank.xmat.1D',np.c_[t,t,np.ones(6)],[1,1,0])
    dy=np.c_[2*t+1,-3*t+4];volume('p6support_rank.nii',dy)
    run([BINARY,'-input','p6support_rank.nii','-matrix','p6support_rank.xmat.1D',
         '-penalty','enet','-lambda',.1,'-l1_ratio',.7,'-support_refit',
         '-prewhite_ar1',.4,'-prefix','p6support_rank','-save_nuisance',
         '-save_fitted','-save_residuals'])
    phase6_check('p6support_rank',phase6_fit(duplicate,dy,[1,1,0],
                 np.arange(6),np.zeros(6),.4,'enet',.1,support=True))
    rank_maps=read('p6support_rank','diagnostics')
    assert np.all(rank_maps[9]==2) and np.all(rank_maps[10]==1)
    # Block CV creates training gaps INSIDE a run. A global pre-whitening of
    # Y would leak held-out samples and disagree with these fold-local fits.
    folds=prediction_folds(good,starts,block=7,gap=1)
    weights=np.array([len(te) for _,te in folds])/len(good)
    for kind in ['banded','support']:
        parameters=vectors if kind=='banded' else np.array([.4,.1,.02])
        strength=vectors.mean(axis=1) if kind=='banded' else parameters
        gridflags=[*band,'-band_penalties','p6.penalties'] if kind=='banded' else [
            '-penalty','enet','-lambdas','.4,.1,.02','-l1_ratio',.7,
            '-support_refit','-standardize','-solver_tol',1e-10]
        losses=[]
        for tr,te in folds:
            row=[]
            for parameter in parameters:
                fit=phase6_fit(x[tr],y[tr],roles,good[tr],runs[tr],.55,
                    'banded' if kind=='banded' else 'enet',parameter,groups,
                    support=kind=='support',standardize=True)
                predicted=x[te][:,roles==1]@fit[0]+x[te][:,roles==0]@fit[1]
                row.append(np.mean((y[te]-predicted)**2,axis=0))
            losses.append(row)
        losses=np.array(losses)
        for rule in ['min','1se']:
            selected,mean=prediction_select(losses,weights,-strength,rule=='1se')
            tag=f'p6_cv_{kind}_{rule}'
            args=[*base,*gridflags,'-prewhite_ar1',.55,'-cv_predict',
                  '-cv_block_size',7,'-cv_gap',1,'-cv_rule',rule,'-save_cv']
            run([BINARY,*args,'-prefix',tag,'-batch_size',2,'-jobs',2])
            dg=read(tag,'diagnostics');offset=7 if kind=='banded' else 12
            assert np.array_equal(dg[offset+3],selected)
            assert_allclose(dg[offset+5:],mean,rtol=6e-5,atol=5e-6)
            for v,c in enumerate(selected):
                fit=phase6_fit(x,y[:,v:v+1],roles,good,runs,.55,
                    'banded' if kind=='banded' else 'enet',parameters[c],groups,
                    support=kind=='support',standardize=True)
                for name,expected in zip(['betas','nuisance','fitted','residuals'],fit):
                    assert_allclose(read(tag,name)[:,v],expected[:,0],rtol=7e-5,atol=6e-6)
            run([BINARY,*args,'-prefix',tag+'_serial','-batch_size',5,'-jobs',1])
            for name in ['betas','nuisance','diagnostics','fitted','residuals']:
                assert np.array_equal(read(tag,name),read(tag+'_serial',name))
        # Global normalized scoring uses the same fold fits; a zero-response
        # voxel has zero baseline and contributes no vote to the global index.
        baseline=np.array([np.mean((y[te]-y[tr].mean(axis=0))**2,axis=0)
                           for tr,te in folds])
        pooled=weights@baseline;defined=pooled>0
        global_loss=(losses[:,:,defined]/pooled[defined]).mean(axis=2)
        global_index,_=prediction_select(global_loss[:,:,None],weights,
                                         -strength,False)
        volume('p6mask.nii',np.ones((1,5)))
        tag='p6_global_'+kind
        run([BINARY,*base,*gridflags,'-prewhite_ar1',.55,'-cv_predict',
             '-cv_block_size',7,'-cv_gap',1,'-cv_global','-mask','p6mask.nii',
             '-prefix',tag,'-batch_size',2,'-jobs',2])
        offset=7 if kind=='banded' else 12
        assert np.all(read(tag,'diagnostics')[offset+3]==global_index[0])

    # HRF selection remains unpenalized, even when the winning design will
    # receive a banded or sparse-support refit. Candidate storage must retain
    # the winner's diagnostics and support, not the last design's results.
    second=x.copy();second[:,roles==1]=rng.normal(size=(len(x),4))
    manifest,designs=candidate_manifest('p6_candidates',[x,second],roles,
                                       good,starts,ids=(11,4))
    selection=np.array([np.mean(phase6_fit(design,y,roles,good,runs,.55,
                              'ridge',1)[3]**2,axis=0) for design in designs])
    selected=np.array([min(range(2),key=lambda c:(selection[c,v],[11,4][c]))
                       for v in range(y.shape[1])])
    for kind in ['banded','support']:
        flags=[*band,'-band_penalties','p6.fixed'] if kind=='banded' else [
            '-penalty','enet','-lambda',.08,'-l1_ratio',.7,'-support_refit',
            '-standardize','-solver_tol',1e-10]
        tag='p6_hrf_'+kind
        run([BINARY,'-input','p6.nii','-candidate_manifest',manifest,
             '-candidate_policy','hrf',*flags,'-prewhite_ar1',.55,
             '-prefix',tag,'-save_fitted','-save_residuals','-save_nuisance',
             '-save_candidates','-batch_size',2,'-jobs',2])
        maps,hdr=candidate_maps(tag)
        assert np.array_equal(maps[0],np.array([11,4])[selected])
        assert_allclose(maps[3:],selection,rtol=8e-5,atol=6e-6)
        for v,c in enumerate(selected):
            expected=phase6_fit(designs[c],y[:,v:v+1],roles,good,runs,.55,
                'banded' if kind=='banded' else 'enet',
                vectors[0] if kind=='banded' else .08,groups,
                support=kind=='support',standardize=True)
            for name,value in zip(['betas','nuisance','fitted','residuals'],expected):
                assert_allclose(read(tag,name)[:,v],value[:,0],rtol=8e-5,atol=6e-6)
            if kind=='support':
                assert read(tag,'diagnostics')[9,v]==expected[4][0]

    # Repeat tuning drops columns as well as rows. Group IDs are deliberately
    # run-specific, so one entire group disappears in each training fold.
    rows=np.delete(np.arange(48),[3,19,37]);rs=rows//16
    trials=np.repeat(np.arange(3),2);conditions=np.tile([0,1],3)
    rx=np.zeros((len(rows),9));rr=np.array([1]*6+[0]*3)
    for r in range(3):
        rx[rs==r,2*r:2*r+2]=rng.normal(size=((rs==r).sum(),2))
        rx[rs==r,6+r]=1
    rx=xmat('p6repeat.xmat.1D',rx,rr,rows,48,[0,16,32])
    ry=(rx@rng.normal(size=(9,3))+rng.normal(size=(len(rows),3))*.3).astype(np.float32).astype(float)
    rf=np.zeros((48,3));rf[rows]=ry;volume('p6repeat.nii',rf)
    np.savetxt('p6repeat.groups',np.r_[trials,[-1]*3],fmt='%d')
    pv=np.array([[.02,.3,.7],[1,.08,.4],[.3,.6,.02]])
    np.savetxt('p6repeat.penalties',pv)
    np.savetxt('p6repeat.trials',np.c_[np.arange(6),np.arange(6)+100,
                                      conditions,trials],fmt='%d')
    reference=phase6_fit(rx,ry,rr,rows,rs,.4,'ridge',1)[0]
    loss=np.zeros((3,3));pairs=0
    for r in range(3):
        tr=rs!=r;columns=np.any(rx[tr]!=0,axis=0);training=np.flatnonzero(trials!=r)
        for c in range(3):
            beta=phase6_fit(rx[tr][:,columns],ry[tr],rr[columns],rows[tr],rs[tr],
                           .4,'banded',pv[c],trials[training])[0]
            for local,global_col in enumerate(training):
                for target in np.flatnonzero((trials==r)&(conditions==conditions[global_col])):
                    loss[c]+=(beta[local]-reference[target])**2
        pairs+=sum(np.sum((trials==r)&(conditions==conditions[j])) for j in training)
    loss/=pairs;selected=np.argmin(loss,axis=0)
    run([BINARY,'-input','p6repeat.nii','-matrix','p6repeat.xmat.1D',
         '-penalty','banded','-band_groups','p6repeat.groups',
         '-band_penalties','p6repeat.penalties','-trial_info','p6repeat.trials',
         '-prewhite_ar1',.4,'-prefix','p6repeat','-save_nuisance','-save_fitted',
         '-save_residuals','-save_cv'])
    dg=read('p6repeat','diagnostics')
    assert_allclose(dg[11:],loss,rtol=8e-5,atol=6e-6)
    assert np.array_equal(dg[10],selected)
    for v,c in enumerate(selected):
        expected=phase6_fit(rx,ry[:,v:v+1],rr,rows,rs,.4,'banded',pv[c],trials)
        for name,value in zip(['betas','nuisance','fitted','residuals'],expected):
            assert_allclose(read('p6repeat',name)[:,v],value[:,0],rtol=8e-5,atol=6e-6)

    # Invalid combinations and tables fail before creating outputs.
    failure_case([*base,*band,'-band_penalties','p6.fixed','-fraction',.5])
    failure_case([*base,'-fraction',.5,'-support_refit'],'support_refit')
    failure_case([*base,'-fraction',.5,'-prewhite_ar1',1],'prewhite_ar1')
    failure_case([*base,'-fraction',.5,'-prewhite_ar1','nan'],'prewhite_ar1')
    for content in ['.1 0\n','.1 nan\n','.1 .2\n.3\n','.1 .2\n.1 .2\n']:
        Path('p6.badpen').write_text(content)
        failure_case([*base,*band,'-band_penalties','p6.badpen'])
    for content in ['0 0 1 0 -1 1\n','0 -1 2 0 -1 1\n','0 -1 0 0 -1 0\n']:
        Path('p6.badgroups').write_text(content)
        failure_case([*base,'-penalty','banded','-band_groups','p6.badgroups',
                      '-band_penalties','p6.fixed'])
    failure_case([*base,*band,'-band_penalties','p6.fixed','-memory_mb',.001])
    print('PASS phase6 banded primal/dual, fixed AR1 boundaries, sparse support '
          'refits, fold-local CV/one-SE, original-unit outputs and scheduling')


with tempfile.TemporaryDirectory(prefix='picard-tests-') as folder:
    os.chdir(folder)
    rng=np.random.default_rng(291126)
    n=24; keep=np.ones(n,int);keep[[0,4,13,21]]=0;good=np.where(keep)[0]
    t=np.arange(n);stim=np.column_stack([np.sin(t*.8),np.cos(t*.37)])
    np.savetxt('a.1D',stim[:,0]);np.savetxt('b.1D',stim[:,1]);np.savetxt('censor.1D',keep,fmt='%d')
    run([DECONVOLVE,'-nodata',n,1.5,'-polort','1','-concat','1D: 0 12',
         '-num_stimts','2','-stim_file','1','a.1D','-stim_label','1','A',
         '-stim_file','2','b.1D','-stim_label','2','B','-censor','censor.1D',
         '-x1D','real.xmat.1D','-x1D_stop'])
    x=np.loadtxt('real.xmat.1D').astype(np.float32).astype(float)
    assert len(x)==len(good) and x.shape[1]==6
    roles=[0,0,0,0,1,1]
    yy=(x@rng.normal(size=(6,5))+rng.normal(size=(len(x),5))*.1).astype(np.float32).astype(float)
    yy[:,3]=0
    full=np.zeros((n,5),np.float32);full[good]=yy;full[keep==0]=np.nan
    volume('input.nii',full);volume('mask.nii',[[1,1,1,1,0]])
    base=['-input','input.nii','-matrix','real.xmat.1D','-mask','mask.nii',
          '-save_nuisance','-save_fitted','-save_residuals','-write_manifest']
    for fraction in (0,.35,1):
        prefix=f'f{fraction}'
        run([BINARY,*base,'-fraction',fraction,'-prefix',prefix,'-batch_size','2','-jobs','2'])
        check_fit(prefix,x,yy,roles,fraction,[1,1,1,1,0])
    rows=np.loadtxt('f0.35.rows.1D',dtype=int)
    assert np.array_equal(rows[:,1],good) and np.array_equal(rows[:,2],good//12)
    manifest=list(csv.DictReader(open('f0.35.betas.tsv'),delimiter='\t'))
    assert [int(r['XColumn']) for r in manifest]==[4,5]
    assert [r['Label'] for r in manifest]==['A#0','B#0']
    labels=run([Path(DECONVOLVE).with_name('3dinfo'),'-label','f0.35+orig.HEAD'])
    assert 'A#0_Coef|B#0_Coef' in labels
    run([BINARY,*base,'-fraction','.35','-prefix','serial','-batch_size','5','-jobs','1'])
    for kind in ('betas','nuisance','diagnostics','fitted','residuals'):
        assert np.array_equal(read('serial',kind),read('f0.35',kind)),kind
    before=Path('f0.35+orig.HEAD').read_bytes()
    run([BINARY,*base,'-fraction','.35','-prefix','f0.35'],ok=False)
    assert Path('f0.35+orig.HEAD').read_bytes()==before
    msg=run([BINARY,'-matrix','real.xmat.1D','-validate_only','-prefix','validate'])
    assert 'Penalized=2 Nuisance=4' in msg and not list(Path('.').glob('validate.*'))
    header=nib.load('f0.35+orig.HEAD').header.info
    assert header['PICARD_ROW_IDS']==good.tolist()
    assert header['PICARD_RUN_IDS']==(good//12).tolist()
    assert header['PICARD_COLUMN_ROLES']==roles
    assert 'minimum-norm' in header['PICARD_MODEL']
    # The default is one native bucket, with all provenance in its header.
    run([BINARY,'-input','input.nii','-matrix','real.xmat.1D','-fraction','.35','-prefix','compact'])
    assert sorted(p.name for p in Path('.').glob('compact*')) in [
        ['compact+orig.BRIK','compact+orig.HEAD'],['compact+orig.BRIK.gz','compact+orig.HEAD']]
    assert len(nib.load('compact+orig.HEAD').header.info['BRICK_LABS'].split('~'))>=9
    # Explicit NIfTI and familiar AFNI fitts/errts prefixes remain available.
    run([BINARY,'-input','input.nii','-matrix','real.xmat.1D','-fraction','.35',
         '-prefix','explicit_nifti.nii.gz','-fitts','my_fitts','-errts','my_errts'])
    assert_allclose(np.asarray(nib.load('explicit_nifti.nii.gz').dataobj).reshape(5,-1),
                    np.asarray(nib.load('compact+orig.HEAD').dataobj).reshape(5,-1),rtol=0,atol=0)
    assert nib.load('my_fitts+orig.HEAD').shape[-1]==len(good)
    assert nib.load('my_errts+orig.HEAD').header.info['PICARD_ROW_IDS']==good.tolist()
    print('PASS real Deconvolve xmat, native bucket, header provenance, censor/run maps, labels, mask, endpoints, batching and workers')

    # Redundant nuisance, nonorthogonal predictors and rank-deficient wide fits.
    for case,(nr,npred,nz) in enumerate([(9,3,2),(5,10,0),(7,4,3),(3,1,0)]):
        z=np.ones((nr,nz));xx=rng.normal(size=(nr,npred))
        if nz>1:z[:,1]=2
        if nz>2:z[:,2]=np.arange(nr)
        if npred>1:xx[:,-1]=xx[:,0]
        matrix=xmat('synthetic.xmat.1D',np.column_stack([z,xx]),[0]*nz+[1]*npred)
        y=rng.normal(size=(nr,4)).astype(np.float32).astype(float);y[:,3]=0
        volume('synthetic.nii',y)
        args=['-input','synthetic.nii','-matrix','synthetic.xmat.1D','-save_fitted','-save_residuals']
        if nz:args+=['-save_nuisance']
        for fraction in (.2,.8,1):
            prefix=f'case{case}_{fraction}'
            run([BINARY,*args,'-fraction',fraction,'-prefix',prefix])
            check_fit(prefix,matrix,y,[0]*nz+[1]*npred,fraction)
    # Effective design entirely in the nuisance span is rank zero.
    matrix=xmat('zero.xmat.1D',np.ones((5,2)),[0,1]);volume('zero.nii',np.arange(10).reshape(5,2))
    run([BINARY,'-input','zero.nii','-matrix','zero.xmat.1D','-fraction','.4','-prefix','zero','-save_nuisance'])
    dg=read('zero','diagnostics');assert np.all(dg[5]==0) and np.all(dg[6]==1) and np.all(read('zero','betas')==0)
    # AFNI brick scaling is applied before the single conversion to double.
    matrix=xmat('scaled.xmat.1D',np.column_stack([np.ones(6),np.arange(6)]),[0,1])
    raw=np.arange(12,dtype=np.int16).reshape(6,2);volume('scaled.nii',raw,np.int16,2.5)
    run([BINARY,'-input','scaled.nii','-matrix','scaled.xmat.1D','-fraction','.6','-prefix','scaled',
         '-save_nuisance','-save_fitted','-save_residuals'])
    check_fit('scaled',matrix,raw*2.5,[0,1],.6)
    # Compare with a separate AFNI program's loaded data, not the original
    # NIfTI bytes. Float64 conversion, byte order, gzip and slope/intercept
    # handling belong to AFNI; Picard must inherit those choices.
    raw64=1e9+np.arange(12,dtype=float).reshape(6,2)*.001
    header=nib.Nifti1Header(endianness='>');header.set_data_dtype(np.float64)
    image=nib.Nifti1Image(raw64.T[:,None,None,:],np.eye(4),header=header)
    image.header.set_slope_inter(1.5,4)
    nib.save(image,'double.nii.gz')
    run([Path(DECONVOLVE).with_name('3dcopy'),'double.nii.gz','double_native'])
    for source, prefix in [('double.nii.gz','double'),
                           ('double_native+orig','double_nativefit')]:
        run([BINARY,'-input',source,'-matrix','scaled.xmat.1D','-fraction','.6',
             '-prefix',prefix,'-save_nuisance','-save_fitted','-save_residuals'])
    for kind in ('betas','nuisance','diagnostics','fitted','residuals'):
        assert np.array_equal(read('double',kind),read('double_nativefit',kind))
    # Native AFNI input continues to use its unchanged dataset reader.
    run([Path(DECONVOLVE).with_name('3dcopy'),'scaled.nii','native'])
    run([BINARY,'-input','native+orig','-matrix','scaled.xmat.1D','-fraction','.6','-prefix','nativefit',
         '-save_nuisance','-save_fitted','-save_residuals'])
    for kind in ('betas','nuisance','diagnostics','fitted','residuals'):
        assert np.array_equal(read('nativefit',kind),read('scaled',kind))
    print('PASS rank deficiency, wide matrices, redundant nuisance, zero reference and scaled input')

    # Nonfinite samples and mask values follow the standard AFNI loader.
    # Use 3dcopy as the I/O reference and NumPy as the numerical reference.
    # Nonzero intercept scaling also exercises the NIfTI conversion path.
    dirty=np.arange(18,dtype=np.float32).reshape(6,3)
    dirty[1,0]=np.nan; dirty[2,1]=np.inf; dirty[3,2]=-np.inf
    image=nib.Nifti1Image(dirty.T[:,None,None,:],np.eye(4))
    image.header.set_slope_inter(1.5,4)
    nib.save(image,'repair.nii.gz')
    volume('repair_mask.nii',[[1,np.nan,1]])
    copy=Path(DECONVOLVE).with_name('3dcopy')
    run([copy,'repair.nii.gz','repair_native'])
    run([copy,'repair_mask.nii','repair_mask_native'])
    loaded=np.asarray(nib.load('repair_native+orig.HEAD').dataobj).reshape(3,6).T
    loaded_mask=np.asarray(nib.load('repair_mask_native+orig.HEAD').dataobj).ravel()!=0
    for source, mask, prefix in [
            ('repair.nii.gz','repair_mask.nii','repair'),
            ('repair_native+orig','repair_mask_native+orig','repair_nativefit')]:
        run([BINARY,'-input',source,'-mask',mask,'-matrix','scaled.xmat.1D',
             '-fraction','.6','-prefix',prefix,'-save_nuisance',
             '-save_fitted','-save_residuals'])
        check_fit(prefix,matrix,loaded,[0,1],.6,loaded_mask)
    for kind in ('betas','nuisance','diagnostics','fitted','residuals'):
        assert np.array_equal(read('repair',kind),read('repair_nativefit',kind))

    # A native float brick with NaN distinguishes the environment policies.
    # Build the header with AFNI, then inject NaN into this disposable fixture.
    # With scanning disabled, Picard must reject the value that reaches it.
    volume('scan.nii',np.arange(12).reshape(6,2))
    run([copy,'scan.nii','scan_native'])
    header=nib.load('scan_native+orig.HEAD').header
    dtype=header.get_data_dtype()
    assert dtype.kind=='f' and dtype.itemsize==4
    brick=Path('scan_native+orig.BRIK')
    compressed=Path(str(brick)+'.gz')
    if compressed.exists():
        brick.write_bytes(gzip.decompress(compressed.read_bytes()))
        compressed.unlink()
    payload=np.frombuffer(brick.read_bytes(),dtype=dtype).copy()
    payload[2]=np.nan
    brick.write_bytes(payload.tobytes())
    scan_args=['-input','scan_native+orig','-matrix','scaled.xmat.1D',
               '-fraction','.6']
    run([BINARY,*scan_args,'-prefix','scan_yes'])
    prior=ENV['AFNI_FLOATSCAN']
    try:
        ENV['AFNI_FLOATSCAN']='NO'
        failure_case(scan_args,'nonfinite',prefix='scan_no')
    finally:
        ENV['AFNI_FLOATSCAN']=prior
    print('PASS AFNI loader parity, double conversion, NaN/Inf repair, masks '
          'and AFNI_FLOATSCAN settings')

    # Explicit role override on a matrix without metadata; no implicit intercept.
    matrix=xmat('plain.xmat.1D',np.column_stack([np.ones(6),np.arange(6)]))
    run([BINARY,'-input','scaled.nii','-matrix','plain.xmat.1D','-penalty_cols','1','-fraction','.6','-prefix','explicit',
         '-save_nuisance','-save_fitted','-save_residuals'])
    for kind in ('betas','nuisance','diagnostics','fitted','residuals'):
        assert np.array_equal(read('explicit',kind),read('scaled',kind))
    default=['-input','scaled.nii','-matrix','plain.xmat.1D','-fraction','.6']
    failure_case(default,'ColumnGroups')
    failure_case(default+['-penalty_cols','0,0'],'duplicate')
    failure_case(default+['-penalty_cols','2'],'out-of-range')
    failure_case(default+['-penalty_cols','1','-nuisance_cols','1'],'disjoint')
    failure_case(default+['-penalty_cols','1','-nuisance_cols','none'],'cover')
    failure_case(base+['-fraction','nan'],'fraction')
    failure_case(base+['-fraction','1.1'],'fraction')
    failure_case(base+['-fraction','.5','-batch_size','0'],'batch_size')
    failure_case(base+['-fraction','.5','-rank_tol','0'],'rank_tol')
    failure_case(base+['-fraction','.5','-memory_mb','.1'],'memory_mb')
    failure_case(['-input','scaled.nii','-matrix','real.xmat.1D','-fraction','.5'],'NRowFull')
    xmat('ambiguous.xmat.1D',np.ones((4,1)),[1],nfull=6)
    failure_case(['-input','scaled.nii','-matrix','ambiguous.xmat.1D','-fraction','.5'],'ambiguous')
    xmat('duplicate.xmat.1D',np.ones((6,2)),[0,1],labels=['same','same'])
    failure_case(['-input','scaled.nii','-matrix','duplicate.xmat.1D','-fraction','.5'],'duplicate')
    volume('empty.nii',np.zeros((1,5)))
    failure_case(['-input','input.nii','-matrix','real.xmat.1D','-fraction','.5','-mask','empty.nii'],'empty')
    img=nib.load('mask.nii');shift=img.affine.copy();shift[0,3]=5
    nib.save(nib.Nifti1Image(np.asarray(img.dataobj),shift),'shift.nii')
    failure_case(['-input','input.nii','-matrix','real.xmat.1D','-fraction','.5','-mask','shift.nii'],'grid')
    # Metadata parsing must not silently fall back to an uncensored timeline.
    original=Path('scaled.xmat.1D').read_text()
    for key,value in [('GoodList','nonsense'),('RunStart','bad'),('ColumnGroups','0,1.5')]:
        text=original
        if key=='ColumnGroups':text=text.replace('ColumnGroups = "0,1"',f'ColumnGroups = "{value}"')
        else:text=text.replace('# >',f'# {key} = "{value}"\n# >')
        if key=='GoodList':text=text.replace('# >','# NRowFull = "6"\n# >')
        Path('malformed.xmat.1D').write_text(text)
        failure_case(['-input','scaled.nii','-matrix','malformed.xmat.1D','-fraction','.5'])
    Path('collision.rows.1D').write_text('keep me')
    run([BINARY,*base,'-fraction','.5','-prefix','collision'],ok=False)
    assert Path('collision.rows.1D').read_text()=='keep me' and not Path('collision+orig.HEAD').exists()
    failure_case(base+['-fraction','.5'],prefix='missing_directory/out')
    # A later companion write failure rolls back the already-written bucket,
    # while preserving the unrelated file that obstructed the output path.
    Path('blocked').write_text('must remain')
    failure_case(['-input','scaled.nii','-matrix','scaled.xmat.1D','-fraction','.5',
                  '-fitts','blocked/fitts'],prefix='rollback')
    assert Path('blocked').read_text()=='must remain'
    print('PASS role validation, input failures, output protection and memory budget')

    if FIXTURE.exists():
        ref=json.loads(FIXTURE.read_text())
        matrix=xmat('fracridge.xmat.1D',ref['x'],[1]*len(ref['x'][0]))
        y=np.asarray(ref['y'],np.float32).astype(float);volume('fracridge.nii',y)
        for i,fraction in enumerate(ref['fractions']):
            prefix=f'reference{i}'
            run([BINARY,'-input','fracridge.nii','-matrix','fracridge.xmat.1D','-fraction',fraction,
                 '-prefix',prefix,'-save_fitted','-save_residuals'])
            check_fit(prefix,matrix,y,[1]*matrix.shape[1],fraction)
            # Fracridge interpolates a log-alpha grid; Picard solves the norm
            # equation directly. A 2% reference tolerance is deliberate.
            assert_allclose(read(prefix,'betas'),np.asarray(ref['coef'])[:,i,:],rtol=.02,atol=.003)
            assert_allclose(read(prefix,'diagnostics')[2],np.asarray(ref['alpha'])[i],rtol=.02,atol=.003)
        print('PASS pinned nrdg/fracridge comparison:',ref['commit'])
    else:
        raise AssertionError('missing pinned fracridge reference fixture')
    repeated_condition_tests(rng)
    prediction_tests(rng)
    elastic_tests(rng)
    candidate_tests(rng)
    phase6_tests(rng)
print('PASS 3dPicard fixed/tuned integration and independent numerical oracles')
