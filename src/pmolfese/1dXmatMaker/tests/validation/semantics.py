# Developed at the NIH
#
# Part of AFNI. See the repository LICENSE.txt for licensing terms and
# disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
# Third-party and adapted portions retain their existing attribution and license
# terms; see THIRD_PARTY.md.

"""Timing, modulation, contrasts, censoring and publication contracts.

Use tests/run.py for the public commands. Fixtures are grouped by reference
implementation and behavior; checks never rewrite expected numerical results.
"""

import argparse
import json
import math
import tempfile
import shutil
import os
import subprocess
from . import ROOT
from pathlib import Path
from .reference import digest, read_matrix, run_case, integers
from .designs import invoke, compare_events as compare, afni_run

R = ROOT


# Full-design comparisons (M4)
# ------------------------------------------------------------------------
CASES=R/'fixtures/afni/semantics/cases.json'
REF=R/'fixtures/afni/semantics/reference'

def floats(text):
    """Decode AFNI's repeated-value notation for numeric contrast attributes."""
    out=[]
    for token in text.split(','):
        if '@' in token:
            n,v=token.split('@');out.extend([float(v)]*int(n))
        else:out.append(float(token))
    return out

def check(actual,expected):
    """Use existing M2 tolerances; normalize only equivalent GLT text encoding."""
    a,_=read_matrix(actual);e,_=read_matrix(expected)
    text=actual.read_text()
    for key in e:
        if key.startswith('GltMatrix_'):
            av,ev=floats(a[key]),floats(e[key])
            assert len(av)==len(ev) and av[:2]==ev[:2],(key,av,ev)
            # AFNI serializes GLT weights using six significant digits. Keep
            # zero locations exact; only nonzero printing precision is tolerant.
            assert all((v==0 if w==0 else math.isclose(v,w,rel_tol=5e-6,abs_tol=0)) for v,w in zip(av[2:],ev[2:])),(key,av,ev)
            text=text.replace(f'{key} = "{a[key]}"',f'{key} = "{e[key]}"')
    # compare() reads a temporary copy solely to ignore compression syntax.
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/'X';p.write_text(text);return compare(p,expected)

def semantics_main(argv=None):
    """Compare full timing/modulation designs, or explicitly capture a new corpus."""
    global CASES,REF
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--capture-reference',type=Path);p.add_argument('--maker',type=Path);p.add_argument('--reader',type=Path)
    p.add_argument('--report',type=Path);p.add_argument('--cases',type=Path,default=CASES);p.add_argument('--reference-dir',type=Path,default=REF)
    a=p.parse_args(argv);CASES=a.cases;REF=a.reference_dir;cases=json.loads(CASES.read_text())
    if a.capture_reference:
        exe=a.capture_reference.resolve(strict=True)
        pinned=json.loads((R/'fixtures/afni/compatibility/reference/manifest.json').read_text())
        assert digest(exe)==pinned['executable_sha256']
        assert not REF.exists(),'Never overwrite reference captures'
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);records=[run_case(exe,c,d/c['id']) for c in cases]
            (d/'manifest.json').write_text(json.dumps(dict(executable_sha256=digest(exe),cases_sha256=digest(CASES),cases=records),indent=2)+'\n')
            shutil.copytree(d,REF)
        return
    manifest=json.loads((REF/'manifest.json').read_text());assert digest(CASES)==manifest['cases_sha256']
    for c in manifest['cases']:
        for name,sha in c['files'].items():assert digest(REF/c['id']/name)==sha
    exe=a.maker.resolve(strict=True);results=[]
    with tempfile.TemporaryDirectory() as d:
        root=Path(d)
        for c in cases:
            work=root/c['id'];work.mkdir()
            for name,text in c['inputs'].items():(work/name).write_text(text)
            invoke(exe,c['maker_args']+['-x1D','X.xmat.1D','-x1D_uncensored','X.full.xmat.1D'],work)
            worst=max(check(work/f,REF/c['id']/f) for f in ['X.xmat.1D','X.full.xmat.1D'])
            if c['id']=='regcensored':
                at,full=read_matrix(work/'X.full.xmat.1D');rt,reg=read_matrix(work/'R.xmat.1D')
                _,refreg=read_matrix(REF/c['id']/'R.xmat.1D')
                good=integers(read_matrix(work/'X.xmat.1D')[0]['GoodList']);bad=[i for i in range(len(full)) if i not in good]
                assert integers(rt['ColumnGroups'])==integers(at['ColumnGroups'])+[0]*len(bad)
                assert len(set(x.strip() for x in rt['ColumnLabels'].split(';')))==len(reg[0])
                for i,row in enumerate(reg):
                    assert row==full[i]+[float(i==j) for j in bad]
                    assert all(math.isclose(v,w,abs_tol=2e-6,rel_tol=2e-5) for v,w in zip(row,refreg[i]))
            if a.reader:
                reader=a.reader.resolve(strict=True)
                afni_run([reader,'-infile',work/'X.xmat.1D','-uncensored',work/'X.full.xmat.1D','-describe'],work)
                if c['id']=='regcensored':afni_run([reader,'-infile',work/'R.xmat.1D','-describe'],work)
            results.append(dict(id=c['id'],max_abs_error=worst,afni_reader=bool(a.reader)))
            print(c['id'],'OK')
    report=dict(maker_sha256=digest(exe),reference_sha256=manifest['executable_sha256'],cases=results)
    if a.report:a.report.write_text(json.dumps(report,indent=2)+'\n')
    print(f'{len(results)} full-design comparisons passed')


# Failure, environment and original-fixture checks (M4)
# ------------------------------------------------------------------------

def check_contracts(exe,root):
    base=['-run_lengths','12','-TR','1','-polort','0','-num_stimts','1','-local_times']
    source=['-stim_times_AM2','1','times','GAM','-stim_label','1','event']
    outputs=['-x1D','X','-x1D_uncensored','F','-x1D_regcensored','R','-event_manifest','events']
    work=root/'boundaries';work.mkdir();count=0
    def fail(args,text='1*1 6*3\n'):
        nonlocal count
        (work/'times').write_text(text)
        invoke(exe,args+outputs,work,2)
        assert not any((work/p).exists() for p in ['X','F','R','events'])
        assert not list(work.glob('.1dXmatMaker-*.tmp'));count+=1
    for text in ['1*1,2 6*3','1*NaN','1*inf','1*1:2','-2*1','12*1','1*1 *','1*1\n6*3\n']:
        fail(base+source,text)
    for options in [['-basis_normall','0'],['-basis_normall','NaN'],['-stim_times_subtract','inf'],['-stim_times_subtract','10']]:
        fail(base+options+source)
    for suffix in [['-CENSORTR','1:12'],['-CENSORTR','2:0'],['-nfirst','12'],['-nlast','12'],['-num_glt','2'],['-gltsym','SYM: missing'],['-gltsym','SYM: event[99]'],['-glt_label','1','missing'],['-stim_nptr','1','2'],['-noblock']]:
        fail(base+source+suffix)
    for center in [':1:2',':NaN',':bad']:
        fail(base+source[:4]+[center]+source[4:])
    fail(base+['-rawAM2']+source[:4]+[':1']+source[4:])
    fail(base+['-allzero_OK','-stim_times_IM','1','times','GAM','-stim_label','1','event'],'*\n')
    for text in ['1 2 3 4','1 -2 3','1 NaN 3']:
        fail(base+['-stim_times_FSL','1','times','dmBLOCK','-stim_label','1','event'],text)
    fail(base+['-stim_times_FSL','1','times','GAM','-stim_label','1','event'],'1 2 3')
    for suffix in [['-stim_nptr','1','0'],['-stim_minlag','1','3','-stim_maxlag','1','1']]:
        fail(base+['-stim_file','1','1D: 12@1','-stim_label','1','s']+suffix)
    # Even the last optional artifact participates in preflight/rollback.
    (work/'times').write_text('1*1 6*3\n');(work/'occupied').write_text('preserve')
    invoke(exe,base+source+['-x1D','X','-x1D_uncensored','F','-x1D_regcensored','occupied','-event_manifest','events'],work,2)
    assert (work/'occupied').read_text()=='preserve' and not any((work/p).exists() for p in ['X','F','events']);count+=1
    # Environmental rawAM2/extend/nodup settings must not alter scientific output.
    matrices=[]
    for index,settings in enumerate([{}, {'AFNI_3dDeconvolve_rawAM2':'YES','AFNI_3dDeconvolve_extend':'YES','AFNI_3dDeconvolve_nodup':'YES'}]):
        env=dict(os.environ,**settings)
        result=subprocess.run([str(exe)]+base+source+['-x1D',f'env{index}'],cwd=work,env=env,capture_output=True,text=True,timeout=30)
        assert result.returncode==0,result.stderr
        matrices.append(read_matrix(work/f'env{index}')[1])
    assert matrices[0]==matrices[1];count+=1
    # Reuse original independent M0 bytes, with no recapture or numerical edits.
    reused=[]
    for c in json.loads((R/'fixtures/afni/compatibility/cases.json').read_text()):
        if not c['id'].startswith(('am','dm','fsl','fs1','local_two','global_two','millisecond','offset','censortr')) or c['id']=='fsl_substitution':continue
        args=c['args'].copy()
        if c['id']=='am2_raw':args.insert(0,'-rawAM2')
        d=root/c['id'];d.mkdir()
        for name,text in c['inputs'].items():(d/name).write_text(text)
        invoke(exe,args+['-x1D','X','-x1D_uncensored','F'],d)
        for actual,expected in [('X','X.xmat.1D'),('F','X.full.xmat.1D')]:check(d/actual,R/'fixtures/afni/compatibility/reference'/c['id']/expected)
        reused.append(c['id'])
    return dict(boundary_and_environment_checks=count,m0_reused=reused)

def contracts_main(argv=None):
    """Check input rejection, publication rollback and environmental independence."""
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--maker',type=Path,required=True);p.add_argument('--report',type=Path);a=p.parse_args(argv)
    with tempfile.TemporaryDirectory() as d:report=check_contracts(a.maker.resolve(strict=True),Path(d))
    if a.report:a.report.write_text(json.dumps(report,indent=2)+'\n')
    print(f"PASS: {report['boundary_and_environment_checks']} contract checks, {len(report['m0_reused'])} original M0 matrix pairs")

