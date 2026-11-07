# Developed at the NIH
#
# Part of AFNI. See the repository LICENSE.txt for licensing terms and
# disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
# Third-party and adapted portions retain their existing attribution and license
# terms; see THIRD_PARTY.md.

"""Pinned GLMsingle and custom sampled-HRF validation and explicit capture.

Use tests/run.py for the public commands. Fixtures are grouped by reference
implementation and behavior; checks never rewrite expected numerical results.
"""

import argparse
import csv
import hashlib
import io
import json
import subprocess
import tempfile
import importlib.util
import shutil
import sys
from . import ROOT, FIXTURES
from pathlib import Path
from .designs import invoke, afni_run
from .reference import read_matrix, integers, digest

R = ROOT


# Read-only waveform/design validation (M6)
# ------------------------------------------------------------------------

def events(tr,duration):
    rows=[(100+i,i%2,i//3,(2+20*(i%3)+.25),duration,f'item{i%2}') for i in range(6)]
    text='trial\tcondition\trun\tonset\tduration\titem\n'+'\n'.join('\t'.join(map(str,row)) for row in rows)+'\n'
    return text,rows

def sampled_main(argv=None):
    """Check pinned curves, full designs, identities and custom-duration failures."""
    import numpy as np
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--maker',type=Path,required=True);p.add_argument('--probe',type=Path,required=True);p.add_argument('--picard',type=Path);p.add_argument('--report',type=Path);a=p.parse_args(argv)
    exe=a.maker.resolve(strict=True);probe=a.probe.resolve(strict=True)
    ref=json.loads((FIXTURES / 'glmsingle/reference.json').read_text());records=[];curves=0;matrices=0;boundaries=0
    provenance=json.loads((R/'data/glmsingle/provenance.json').read_text())
    assert ref['provenance']==provenance
    for name,expected_hash in provenance['source_sha256'].items():
        directory=R/'data/glmsingle' if name=='getcanonicalhrflibrary.tsv' else FIXTURES / 'glmsingle/upstream'
        assert digest(directory/name)==expected_hash,name
    with tempfile.TemporaryDirectory(prefix='maker_m6_') as d:
        root=Path(d)
        for k,case in enumerate(ref['cases']):
            duration,tr=case['duration'],case['tr']
            for mode in ['peak','library']:
                result=subprocess.run([str(probe),str(duration),str(tr),mode],capture_output=True,text=True,check=True)
                actual=np.loadtxt(io.StringIO(result.stdout)).T;expected=np.asarray(case[mode])
                np.testing.assert_array_equal(actual,expected);curves+=20
            work=root/f'case{k}';work.mkdir();text,trials=events(tr,duration);(work/'events').write_text(text)
            nrun=max(80,int(np.ceil(80/tr)));nfull=2*nrun
            keep=np.ones(nfull,dtype=int);keep[[0,3,nrun-1,nrun,nrun+23,nfull-1]]=0;np.savetxt(work/'keep',keep,fmt='%d')
            args=['-candidate_library','glmsingle','-events','events','-run_lengths',f'{nrun},{nrun}','-TR',str(tr),'-polort','0','-censor','keep','-candidate_prefix','set']
            invoke(exe,args,work)
            catalog=list(csv.DictReader((work/'set.candidates.tsv').read_text().splitlines()[1:],delimiter='\t'))
            assert all(float(row['sampling'])==.1 for row in catalog)
            labels=list(csv.DictReader((work/'set.trials.columns.tsv').open(),delimiter='\t'))
            assert [int(r['trial']) for r in labels if r['role']=='penalized']==list(range(100,106))
            for candidate in range(20):
                hrf=np.asarray(case['peak'][candidate]);prepared=np.loadtxt(work/f'set.id{candidate}.hrf.1D')
                np.testing.assert_array_equal(prepared,hrf)
                attrs,retained=read_matrix(work/f'set.id{candidate}.trials.xmat.1D');_,full=read_matrix(work/f'set.id{candidate}.trials.full.xmat.1D');full=np.asarray(full)
                expected=np.zeros((nfull,8));expected[:nrun,0]=1;expected[nrun:,1]=1
                for column,(_,_,run,onset,_,_) in enumerate(trials,2):
                    times=np.arange(nrun)*tr-onset
                    expected[run*nrun:(run+1)*nrun,column]=np.interp(times,np.arange(len(hrf))*tr,hrf,left=0,right=0)
                np.testing.assert_allclose(full,expected,atol=3e-6,rtol=2e-5)
                np.testing.assert_array_equal(retained,full[keep==1]);assert integers(attrs['GoodList'])==np.flatnonzero(keep).tolist()
                _,onoff=read_matrix(work/f'set.id{candidate}.onoff.full.xmat.1D');onoff=np.asarray(onoff)
                np.testing.assert_allclose(onoff[:,2],full[:,2:].sum(axis=1),atol=3e-6,rtol=2e-5)
                np.testing.assert_array_equal(onoff[:,:2],full[:,:2]);matrices+=2
            if a.picard:
                for catalog in ['set.candidates.tsv','set.onoff.candidates.tsv']:
                    afni_run([a.picard.resolve(),'-candidate_manifest',catalog,'-candidate_policy','hrf','-validate_only'],work)
            # Byte-level repeated generation with identical argv in another cwd.
            repeat=root/f'repeat{k}';repeat.mkdir();(repeat/'events').write_text(text);np.savetxt(repeat/'keep',keep,fmt='%d');invoke(exe,args,repeat)
            for path in work.glob('set*'):assert path.read_bytes()==(repeat/path.name).read_bytes()
            records.append(dict(duration=duration,tr=tr,candidates=20,repeat_bytes='identical'))
        work=root/'custom';work.mkdir();text,_=events(1,1);(work/'events').write_text(text)
        (work/'impulse').write_text('0\n1\n0\n');(work/'convolved').write_text('0\n.5\n.5\n0\n')
        header='id\twaveform\tdt\tkind\tduration\tsource\tversion\n'
        library=header+'7\timpulse\t1\timpulse\t0\tanalytic-triangle\tv1\n3\tconvolved\t1\tconvolved\t1\tanalytic-box-integral\tv1\n'
        (work/'library').write_text(library)
        args=['-candidate_library','library','-events','events','-run_lengths','80,80','-TR','1','-polort','0','-hrf_normalization','none','-candidate_prefix','custom']
        invoke(exe,args,work)
        assert read_matrix(work/'custom.id7.trials.xmat.1D')[1]==read_matrix(work/'custom.id3.trials.xmat.1D')[1]
        np.testing.assert_array_equal(np.loadtxt(work/'custom.id7.hrf.1D'),[0,.5,.5,0])
        # Bad data must fail before ANY member of the collection is published.
        for invalid in [library.replace('\tconvolved\t1','\tconvolved\t2'),library.replace('7\t','3\t'),library.replace('\t1\timpulse','\t0\timpulse'),library.replace('\timpulse\t0','\tunknown\t0')]:
            (work/'bad_library').write_text(invalid);bad=args.copy();bad[1]='bad_library';bad[-1]='bad';invoke(exe,bad,work,2);assert not list(work.glob('bad.*'));boundaries+=1
        for option in [['-candidate_ids','7'],['-candidate_ids','7,999'],['-candidate_ids','7,7']]:
            bad=args[:-1]+['bad']+option;invoke(exe,bad,work,2);assert not list(work.glob('bad.*'));boundaries+=1
        # An occupied late catalog prevents publication of earlier matrices.
        (work/'occupied.candidates.tsv').write_text('preserve');invoke(exe,args[:-1]+['occupied'],work,2)
        assert list(work.glob('occupied*'))==[work/'occupied.candidates.tsv'];boundaries+=1
    report=dict(result='pass',maker_sha256=digest(exe),probe_sha256=digest(probe),reference_sha256=digest(FIXTURES / 'glmsingle/reference.json'),waveform_curves=curves,full_designs=matrices,cases=records,custom_no_double_convolution='pass',boundary_checks=boundaries)
    if a.report:a.report.write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {curves} upstream curves, {matrices} trial/ONOFF designs, custom convolution and {boundaries} failures')


# Explicit upstream capture; never invoked by the suite
# ------------------------------------------------------------------------
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def capture_main(argv=None):
    """Explicitly capture pinned upstream curves; --help never writes fixtures."""
    parser=argparse.ArgumentParser(description=capture_main.__doc__)
    parser.add_argument('--output',type=Path,required=True,help='new reference JSON outside the immutable fixture cache')
    args=parser.parse_args(argv)
    output=args.output.resolve()
    if output==FIXTURES or FIXTURES in output.parents:
        parser.error('capture output must be outside the immutable fixture cache')
    if output.exists():raise RuntimeError('Reference already exists; use an explicit new version')
    import numpy as np
    provenance=json.loads((R/'data/glmsingle/provenance.json').read_text())
    upstream=FIXTURES / 'glmsingle/upstream'
    for name,digest in provenance['source_sha256'].items():
        path=R/'data/glmsingle'/name if name=='getcanonicalhrflibrary.tsv' else upstream/name
        assert sha(path)==digest,name
    with tempfile.TemporaryDirectory() as d:
        d=Path(d)
        for name in ['gethrf.py','alt_round.py']:shutil.copy2(upstream/name,d/name)
        shutil.copy2(R/'data/glmsingle/getcanonicalhrflibrary.tsv',d/'getcanonicalhrflibrary.tsv')
        def module(name,path):
            spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
        module('glmsingle.utils.alt_round',d/'alt_round.py');oracle=module('glmsingle_reference',d/'gethrf.py')
        cases=[]
        for duration,tr in [(0,.1),(.05,1),(.15,.7),(1.25,2),(2.5,.8),(4,1.3)]:
            raw=oracle.getcanonicalhrflibrary(duration,tr)
            # This is the separate fit() normalization step on float32 data.
            peak=(raw/raw.max(axis=1,keepdims=True)).astype(np.float32)
            cases.append(dict(duration=duration,tr=tr,library=raw.tolist(),peak=peak.tolist()))
        output.write_text(json.dumps(dict(provenance=provenance,cases=cases),separators=(',',':'))+'\n')
    print('Captured six real-library duration/TR grids, all twenty candidates')

