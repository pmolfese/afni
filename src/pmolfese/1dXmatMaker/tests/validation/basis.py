# Developed at the NIH
#
# Part of AFNI. See the repository LICENSE.txt for licensing terms and
# disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
# Third-party and adapted portions retain their existing attribution and license
# terms; see THIRD_PARTY.md.

"""Response catalog, EXPR grammar and independent numerical anchors.

Use tests/run.py for the public commands. Fixtures are grouped by reference
implementation and behavior; checks never rewrite expected numerical results.
"""

import argparse
import json
import subprocess
import tempfile
from . import ROOT
from pathlib import Path
from .reference import digest, read_matrix, run_case
from .designs import compare_events as compare, afni_run, invoke



# Frozen waveform/EXPR comparisons (M3)
# ------------------------------------------------------------------------


CASES = ROOT / 'fixtures/afni/basis/cases.json'
REFERENCE = ROOT / 'fixtures/afni/basis/reference'


def m0_cases():
    """Reuse original oracle bytes; map only already-understood CLI acquisition
    spelling, without editing formulas, labels, timing files or matrix values.
    """
    for case in json.loads((ROOT/'fixtures/afni/compatibility/cases.json').read_text()):
        a=case['args']
        if case['outcome']!='success' or '-stim_times' not in a:
            continue
        if not case['id'].startswith(('gam','twogam','block','igfun','ublock','tent',
                'csplin','poly','sin','trig','spmg','wav','mion','expr')):
            continue
        if case['id']=='spmg_unknown_fallback':
            continue  # Unknown names deliberately remain errors in Maker.
        n,tr=a[a.index('-nodata')+1:a.index('-nodata')+3]
        i=a.index('-stim_times')
        args=['-run_lengths',n,'-TR',tr,'-polort',a[a.index('-polort')+1],
              '-num_stimts','1','-local_times',*a[i:i+4],'-stim_label','1','event']
        yield dict(case, maker_args=args, kind='matrix', duration=None)


def check_expression_contract(probe):
    """Audit all 123 frozen names/argument counts, independent of Rust dispatch.
    Count=0 parses without evaluating a distribution at invented parameters.
    Disabled source routines must fail explicitly; unknown names never fall back.
    """
    arities=json.loads((ROOT/'fixtures/afni/basis/expression_arities.json').read_text())['arities']
    policies=json.loads((ROOT/'fixtures/afni/basis/expression_coverage.json').read_text())
    inventory=json.loads((ROOT/'tests/validation/source_inventory.json').read_text())
    assert set(arities)==set(policies)=={f['name'] for f in inventory['expression_functions']}
    checks=0
    for name,count in arities.items():
        valid_count=count if count>0 else (3 if name=='WITHIN' else 2)
        valid=','.join(['t']*valid_count)
        for arguments,expected in [(valid,policies[name]!='reference-disabled'),
                                    (','.join(['t']*(count+1)) if count>0 else '',False)]:
            r=subprocess.run([str(probe),f'EXPR(0,1) {name}({arguments})','0','1','0'],
                             capture_output=True,text=True)
            assert (r.returncode==0)==expected,(name,arguments,r.stderr)
            if policies[name]=='reference-disabled':assert 'unimplemented' in r.stderr
            checks+=1
    return checks


def check_random_cli(maker, work):
    """Censoring must select an already-generated stream, not resample it.
    The command provenance changes, so compare numbers and metadata separately.
    """
    outputs=[]
    for index,(seed,censor) in enumerate([(123,False),(123,False),(124,False),(123,True)]):
        directory=work/f'random_{index}';directory.mkdir()
        (directory/'times').write_text('0 8\n')
        (directory/'keep').write_text(''.join('0\n' if i%3==0 else '1\n' for i in range(24)))
        arguments=['-run_lengths','24','-TR','1','-polort','-1','-num_stimts','1',
                   '-local_times','-expr_seed',str(seed),'-stim_times','1','times',
                   'EXPR(0,8) uran(t)~gran(t,1)', '-stim_label','1','event',
                   '-x1D','X.xmat.1D','-x1D_uncensored','X.full.xmat.1D']
        if censor:arguments+=['-censor','keep']
        invoke(maker,arguments,directory)
        outputs.append(read_matrix(directory/'X.full.xmat.1D')[1])
        if censor:
            assert read_matrix(directory/'X.xmat.1D')[1]==[r for i,r in enumerate(outputs[-1]) if i%3]
    assert outputs[0]==outputs[1]==outputs[3]
    assert outputs[0]!=outputs[2]
    return 4


def basis_main(argv=None):
    """Check the fixed AFNI waveform catalog, expression grammar and seeded CLI."""
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--capture-reference',type=Path)
    p.add_argument('--reference-dir',type=Path,default=REFERENCE)
    p.add_argument('--maker',type=Path)
    p.add_argument('--probe',type=Path)
    p.add_argument('--reader',type=Path)
    p.add_argument('--report',type=Path)
    a=p.parse_args(argv)
    cases=json.loads(CASES.read_text())
    if a.capture_reference:
        # Requiring the exact M0 executable prevents a changed producer from
        # moving the expected answer to match a bug in the Rust implementation.
        executable=a.capture_reference.resolve(strict=True)
        pin=json.loads((ROOT/'fixtures/afni/compatibility/reference/manifest.json').read_text())
        assert digest(executable)==pin['executable_sha256']
        a.reference_dir.mkdir()  # A capture must always use a new directory.
        records=[]
        for case in cases:
            records.append(run_case(executable,case,a.reference_dir/case['id']))
        manifest=dict(executable_sha256=digest(executable),cases_sha256=digest(CASES),cases=records)
        (a.reference_dir/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        print(f'CAPTURED: {len(cases)} immutable M3 cases')
        return
    if not a.maker or not a.probe:
        p.error('--maker and --probe are required for validation')
    maker=a.maker.resolve(strict=True);probe=a.probe.resolve(strict=True)
    reader=a.reader.resolve(strict=True) if a.reader else None
    manifest=json.loads((a.reference_dir/'manifest.json').read_text())
    assert manifest['cases_sha256']==digest(CASES)
    for record in manifest['cases']:
        for name,sha in record['files'].items():
            assert digest(a.reference_dir/record['id']/name)==sha
    records=[];failures=[]
    expression_checks=check_expression_contract(probe)
    with tempfile.TemporaryDirectory(prefix='xmat_m3_') as temporary:
        work=Path(temporary)
        random_checks=check_random_cli(maker,work)
        selected=[(c,a.reference_dir/c['id']) for c in cases]
        selected += [(c,ROOT/'fixtures/afni/compatibility/reference'/c['id']) for c in m0_cases()]
        for index,(case,reference) in enumerate(selected):
            directory=work/f'{index}_{case["id"]}';directory.mkdir()
            try:
                if case['id']=='expr_absextreme':
                    # The pinned scalar interpreter truncates its opcode to
                    # eight bytes, so this ten-character dispatch never runs.
                    # Its apparent result is an argument count, not a maximum.
                    # Maker must reject this explicitly instead of emulating
                    # a meaningless value or silently supplying a new model.
                    r=subprocess.run([str(probe), 'EXPR(0,16) absextreme(t,-3,2)',
                                      '0','1','1'],capture_output=True,text=True)
                    assert r.returncode != 0 and 'unimplemented' in r.stderr
                    worst=None
                elif case['kind']=='duration_api':
                    # M3 tests dm waveform math via the public library API.
                    # Married timing/AM matrix assembly remains the M4 gate.
                    args=case['args'];i=args.index('-stim_times_AM1')
                    n=int(args[args.index('-nodata')+1]);tr=args[args.index('-nodata')+2]
                    r=subprocess.run([str(probe),args[i+3],'0',tr,str(n),str(case['duration'])],
                                     capture_output=True,text=True,check=True)
                    actual=[[float(x) for x in line.split()] for line in r.stdout.splitlines()]
                    _,expected=read_matrix(reference/'X.xmat.1D')
                    assert len(actual)==len(expected)
                    worst=0.0
                    for row,(av,ev) in enumerate(zip(actual,expected)):
                        assert len(av)==len(ev)
                        for v,e in zip(av,ev):
                            assert abs(v-e)<=2e-6+2e-5*abs(e),(row,v,e)
                            worst=max(worst,abs(v-e))
                else:
                    for name,text in case['inputs'].items():
                        (directory/name).write_text(text)
                    invoke(maker,case['maker_args']+['-x1D','X.xmat.1D',
                           '-x1D_uncensored','X.full.xmat.1D'],directory)
                    worst=max(compare(directory/name,reference/name)
                              for name in ('X.xmat.1D','X.full.xmat.1D'))
                    if reader:
                        afni_run([reader,'-infile',directory/'X.xmat.1D','-describe'],directory)
                records.append(dict(case=case['id'],kind=case['kind'],max_abs_difference=worst))
            except Exception as error:
                failures.append(dict(case=case['id'],error=str(error)))
    report=dict(maker=str(maker),maker_sha256=digest(maker),probe_sha256=digest(probe),
                expression_contract_checks=expression_checks,random_cli_checks=random_checks,
                cases=records,failures=failures,result='FAIL' if failures else 'pass')
    if a.report:a.report.write_text(json.dumps(report,indent=2)+'\n')
    print(f'{report["result"]}: {len(records)} passing reference cases (including explicit rejection policy); {len(failures)} failures')
    for failure in failures:print(failure)
    if failures:raise SystemExit(1)


# Independent SciPy anchors (M3)
# ------------------------------------------------------------------------
def math_main(argv=None):
    """Compare against independent SciPy distributions, not captured AFNI values."""
    import argparse
    import hashlib
    import json
    import subprocess
    import numpy as np
    import scipy
    from scipy import special as sp, stats
    from pathlib import Path
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe',type=Path,required=True)
    parser.add_argument('--report',type=Path)
    args=parser.parse_args(argv)
    probe=args.probe.resolve(strict=True)
    records=[]
    failures=[]
    def check(expr,x,expected,tol=2e-6):
     formula='EXPR(-100,100) '+expr
     r=subprocess.run([str(probe),formula,str(x), '0','1'],capture_output=True,text=True)
     if r.returncode: failures.append(dict(expression=expr,x=x,error=r.stderr)); return
     v=float(r.stdout);e=float(expected)
     record=dict(expression=expr,x=x,actual=v,expected=e)
     if not np.isfinite(v) or abs(v-e)>tol+2e-6*abs(e): failures.append(record)
     else: records.append(record)
    for n,f in [('j0',sp.j0),('j1',sp.j1),('y0',sp.y0),('y1',sp.y1),('gamma',sp.gamma),('erf',sp.erf),('erfc',sp.erfc),('qg',stats.norm.sf)]:
     for x in [.0001,.1,.5,1,2,5,8,11.9,12,15,20,35]:check(n+'(t)',x,f(x))
    params={2:(10,0,0),3:(10,0,0),4:(3,12,0),5:(0,0,0),6:(4,0,0),7:(2,4,0),8:(12,.4,0),9:(2,.5,0),10:(3,0,0),11:(1,2,0),12:(3,12,2),13:(4,2,0),14:(1,2,0),15:(1,2,0),16:(-2,5,0),17:(8,2,0),18:(1,2,3),19:(4,0,0),20:(2,3,0),21:(1,2,0)}
    dists={2:stats.beta(5,5,loc=-1,scale=2),3:stats.t(10),4:stats.f(3,12),5:stats.norm(),6:stats.chi2(4),7:stats.beta(2,4),8:stats.binom(12,.4),9:stats.gamma(2,scale=2),10:stats.poisson(3),11:stats.norm(1,2),12:stats.ncf(3,12,2),13:stats.ncx2(4,2),14:stats.logistic(1,2),15:stats.laplace(1,2),16:stats.uniform(-2,7),17:stats.nct(8,2),18:stats.weibull_min(3,loc=1,scale=2),19:stats.chi(4),20:stats.invgauss(2/3,scale=3),21:stats.gumbel_r(1,2)}
    for code,d in dists.items():
     a,b,c=params[code]
     for t in [-10,-2,-.1,0,.01,.5,1,2,10,30,80]:
      if code in [8,10] and t!=int(t):continue
      check(f'stat2cdf(t,{code},{a},{b},{c})',t,d.cdf(t))
     if code in [8,10]:continue # AFNI treats the count as continuous for inversion.
     for p in [1e-6,.001,.1,.4,.5,.9,.999,1-1e-6]:
      check(f'cdf2stat({p}+0*t,{code},{a},{b},{c})',0,d.ppf(p),tol=1e-5)

    report=dict(result='FAIL' if failures else 'pass',scipy=scipy.__version__,
        probe_sha256=hashlib.sha256(probe.read_bytes()).hexdigest(),checks=len(records),failures=failures)
    if args.report:args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if failures:raise SystemExit(1)

