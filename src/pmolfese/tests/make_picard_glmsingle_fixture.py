#!/usr/bin/env python3
"""Regenerate the scoring fixture from downloaded, pinned upstream sources.

Usage: python make_picard_glmsingle_fixture.py /path/to/reference_directory
The directory must contain unmodified calcbadness.py and zerodiv.py from
GLMsingle commit 082851f9db87ff30c49370f27365ac1112f18ff7. No package installation
or network access is needed. The normal test consumes only the saved JSON.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np

reference = Path(sys.argv[1])
expected_hashes = {
    'calcbadness.py': 'a87a3a1d7f6ff286a253c3caef0eb1b08c037781cfcd5b942bc5413f2b373d75',
    'zerodiv.py': '645df415d007063cdaec4be3ad07eb2413f022aca7a5bffb5a5ad711a0aa9232',
}
for filename, expected in expected_hashes.items():
    actual = hashlib.sha256((reference / filename).read_bytes()).hexdigest()
    if actual != expected:
        raise ValueError(f'{filename} does not match the pinned source digest')

def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, reference / filename)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result

module('glmsingle.utils.zerodiv', 'zerodiv.py')
scoring = module('picard_reference_scoring', 'calcbadness.py')
# Equal-norm orthogonal trial columns make beta(f)=f*OLS in every fold and
# the full design. Thus this fixture isolates the shared scoring component;
# general designs are tested separately against an independent fold oracle.
betas = np.array([[3, 0, 2, 1, 4, 2, -1, 7],
                  [1, 1, 2, 1, 2, 1, 2, -3],
                  [4, 3, -1, -2, 5, 0, -4, 2]], dtype=float)
conditions = [[0, 0, 1], [0, 1], [0, 1, 9]]
valid = np.empty(3, dtype=object)
stim = np.empty(3, dtype=object)
for i, indices in enumerate([[0, 1, 2], [3, 4], [5, 6, 7]]):
    valid[i] = np.array(indices)
    stim[i] = np.array(conditions[i])
fractions = [1, .75, .5, .25, 0]
results = [betas * f for f in fractions]
# A scalar ndarray matches the helper's single-held-out-run indexing API.
badness = scoring.calcbadness([np.array(i) for i in range(3)], valid, stim,
                             results, np.ones(3, dtype=int))
fixture = dict(
    repository='https://github.com/cvnlab/GLMsingle',
    commit='082851f9db87ff30c49370f27365ac1112f18ff7',
    source_sha256={name: hashlib.sha256((reference/name).read_bytes()).hexdigest()
                   for name in ['calcbadness.py', 'zerodiv.py']},
    scope='calcbadness ordered-pair scoring; one session; no autoscale',
    note='Upstream SSE is in reference sample-SD units. Picard MSE is in '
         'supplied coefficient units. Convert by pair_count/reference_variance. '
         'Equal singular values make full-data and fold-local fractions agree '
         'for this component fixture; no end-to-end equivalence is claimed.',
    betas=betas.T.tolist(), conditions=sum(conditions, []),
    run_ids=[0,0,0,1,1,2,2,2], fractions=fractions,
    badness=badness.T.tolist(), winners=np.argmin(badness, axis=1).tolist())
Path(__file__).with_name('picard_glmsingle_reference.json').write_text(
    json.dumps(fixture, indent=2)+'\n')
