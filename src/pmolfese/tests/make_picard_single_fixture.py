#!/usr/bin/env python3
"""Regenerate a controlled staged fixture using pinned GLMsingle components.

Usage: python make_picard_single_fixture.py REFERENCE_DIR
Supply unmodified glmsingle.py, calcbadness.py, zerodiv.py, and
select_noise_regressors.py from the commit recorded below. Source digests are
checked before executing upstream code. Normal tests read only the saved JSON.

This is NOT an invocation of GLM_single.fit: the controlled design and masks
are supplied, OLS uses NumPy, and the final fraction is 1. It tests the complete
Picard adapter against a staged reference assembled from upstream PCA, scoring
and stopping code, without claiming automatic-policy/full-package equivalence.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import textwrap
import numpy as np
from sklearn.preprocessing import normalize

reference=Path(sys.argv[1])
hashes={
    'glmsingle.py':'878e0db7571dea9f184d9c476429edb5841a22bb4a43afcc25b51419ae8a25dd',
    'select_noise_regressors.py':'ada4e421d7aa5b3a8fc60deec73dd836e3e164a60306c2446f6d8ee04a04c816',
    'calcbadness.py':'a87a3a1d7f6ff286a253c3caef0eb1b08c037781cfcd5b942bc5413f2b373d75',
    'zerodiv.py':'645df415d007063cdaec4be3ad07eb2413f022aca7a5bffb5a5ad711a0aa9232',
}
for name,digest in hashes.items():
    if hashlib.sha256((reference/name).read_bytes()).hexdigest()!=digest:
        raise ValueError(f'{name}: source digest mismatch')
def module(name,filename):
    spec=importlib.util.spec_from_file_location(name,reference/filename)
    obj=importlib.util.module_from_spec(spec);sys.modules[name]=obj
    spec.loader.exec_module(obj);return obj
module('glmsingle.utils.zerodiv','zerodiv.py')
scoring=module('reference_scoring','calcbadness.py')
selection=module('reference_selection','select_noise_regressors.py')
# Execute the actual upstream PCA block, including its float32 projector
# and output casts. Sign alignment is only for comparing equivalent PCs.
source=(reference/'glmsingle.py').read_text()
start=source.index('                noise_pool = polymatrix[run_i].astype')
end=source.index('pcregressors.append(unitary.astype(np.float32))',start)
end+=len('pcregressors.append(unitary.astype(np.float32))')
block=textwrap.dedent(source[start:end])
rng=np.random.default_rng(171)
nr,nt,p,nv=3,16,6,7;n=nr*nt
rows=np.repeat(np.arange(nr),nt)
x=np.zeros((n,p+nr))
noise=rng.normal(size=n)
for r in range(nr):
    x[rows==r,r*2:r*2+2]=rng.normal(size=(nt,2))
    x[rows==r,p+r]=1
    noise[rows==r]+=(r-1)*3*x[rows==r,r*2]
x=x.astype(np.float32).astype(float)
y=np.outer(noise,rng.normal(size=nv))*2
for v in range(3):y[:,v]+=x[:,:p]@np.tile([4,-5],nr)
y+=rng.normal(size=y.shape)*.05
y=y.astype(np.float32).astype(float)
pcs=[]
for r in range(nr):
    context=dict(np=np,normalize=normalize,noise_pool=y[rows==r,3:],
                 polymatrix=[np.eye(nt)-np.ones((nt,nt))/nt],run_i=0,
                 params={'n_pcs':1},pcregressors=[])
    exec(block,context)
    u=context['pcregressors'][0]
    u*=np.sign(u[np.argmax(abs(u),axis=0),np.arange(2)])
    pcs.append(u)
pcs=np.vstack(pcs)
results=[];designs=[]
for count in range(3):
    extra=np.zeros((n,nr*count))
    for r in range(nr):extra[rows==r,r*count:(r+1)*count]=pcs[rows==r,:count]
    design=np.c_[x,extra];designs.append(design)
    results.append(np.linalg.lstsq(design,y[:,:3],rcond=1e-10)[0][:p].T)
valid=np.empty(nr,dtype=object);stim=np.empty(nr,dtype=object)
for r in range(nr):valid[r]=np.arange(r*2,r*2+2);stim[r]=np.arange(2)
badness=scoring.calcbadness([np.array(r) for r in range(nr)],valid,stim,
                          results,np.ones(nr,dtype=int))
trend=np.median(badness,axis=0)
chosen=selection.select_noise_regressors(-trend,1.05)
assert chosen>0
final=np.linalg.lstsq(designs[chosen],y[:,:3],rcond=1e-10)[0]
fixture=dict(repository='https://github.com/cvnlab/GLMsingle',
    commit='082851f9db87ff30c49370f27365ac1112f18ff7',source_sha256=hashes,
    scope='controlled end-to-end Picard case, pinned upstream PCA/scoring/stop; '
          'one supplied HRF, explicit pool and selection masks, one session, '
          'run intercepts, final fraction=1; NOT GLM_single.fit equivalence',
    matrix=x.tolist(),data=y.tolist(),run_length=nt,
    trial_runs=np.repeat(np.arange(nr),2).tolist(),conditions=np.tile([0,1],nr).tolist(),
    pcs=pcs.tolist(),badness=badness.tolist(),pair_count=12,
    chosen=int(chosen),stage_b=results[0].T.tolist(),
    final_beta=final[:p].tolist(),final_fitted=(designs[chosen]@final).tolist())
Path(__file__).with_name('picard_single_reference.json').write_text(json.dumps(fixture,indent=2)+'\n')
print('Wrote pinned staged reference; chosen PC count:',chosen)
