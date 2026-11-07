# Developed at the NIH
#
# Part of AFNI. See the repository LICENSE.txt for licensing terms and
# disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
# Third-party and adapted portions retain their existing attribution and license
# terms; see THIRD_PARTY.md.

"""Baseline, event identity, matrix comparison and CLI-boundary checks.

Use tests/run.py for the public commands. Fixtures are grouped by reference
implementation and behavior; checks never rewrite expected numerical results.
"""

import argparse
import json
import os
import subprocess
import tempfile
import csv
import io
from . import ROOT, FIXTURES
from pathlib import Path
from .reference import digest, integers, read_matrix, run_case, ENV



# Baseline/nuisance matrix contracts (M1)
# ------------------------------------------------------------------------


BASELINE_CASES = FIXTURES / "afni/baselines/cases.json"
BASELINE_REFERENCE = FIXTURES / "afni/baselines/reference"


def compare_baselines(actual, expected):
    """Compare semantic metadata and all elements, not provenance or list spelling.

    AFNI compresses integer lists and Maker writes explicit ones. These encode
    the same identities. Extra XmatMaker provenance does not change the matrix.
    Baseline P0/P1 uses the M0 1e-7 absolute gate; supplied/higher-order values
    use the predeclared 2e-6+2e-5*abs(reference) float-oracle parity threshold.
    """
    aa, ar = read_matrix(actual)
    ea, er = read_matrix(expected)
    assert len(ar) == len(er) and len(ar[0]) == len(er[0])
    integer_keys = {"GoodList", "RunStart", "ColumnGroups"}
    ignored = {"CommandLine", "XmatMakerVersion", "XmatMakerReference"}
    assert aa.keys() - ignored == ea.keys() - ignored, (aa.keys(), ea.keys())
    for key in ea.keys() - ignored:
        if key in integer_keys:
            assert integers(aa[key]) == integers(ea[key]), key
        elif key == "RowTR":
            assert abs(float(aa[key]) - float(ea[key])) < 1e-14
        elif key == "ni_type":
            # NIML permits both 'double' and '1*double' for a single column.
            assert aa[key].removeprefix("1*") == ea[key].removeprefix("1*")
        else:
            assert aa[key] == ea[key], (key, aa[key], ea[key])
    labels = [s.strip() for s in ea["ColumnLabels"].split(";")]
    max_error = 0.0
    for row, (av, ev) in enumerate(zip(ar, er)):
        for col, (a, e) in enumerate(zip(av, ev)):
            tolerance = 2e-6 + 2e-5 * abs(e)
            if labels[col].startswith("Run#") and labels[col].endswith(("Pol#0", "Pol#1")):
                tolerance = 1e-7
            assert abs(a - e) <= tolerance, (actual, row, col, a, e, tolerance)
            # Structural polynomial zeros outside their run must remain exact.
            if e == 0.0 and labels[col].startswith("Run#"):
                starts = integers(ea["RunStart"])
                run = int(labels[col].split("Pol#")[0][4:]) - 1
                full_row = integers(ea["GoodList"])[row]
                end = starts[run + 1] if run + 1 < len(starts) else int(ea["NRowFull"])
                if not starts[run] <= full_row < end:
                    assert a == 0.0
            max_error = max(max_error, abs(a - e))
    return max_error


def invoke(exe, args, cwd, expected=0):
    """Run this exact Maker binary, ensuring errors are ordinary exit codes."""
    result = subprocess.run([str(exe)] + args, cwd=cwd, capture_output=True,
                            text=True, timeout=30)
    assert result.returncode == expected, (args, result.returncode, result.stderr)
    assert "panicked" not in result.stderr
    return result


def negative_and_boundary_checks(exe, work):
    """Invalid input must not publish outputs or damage existing files."""
    work.mkdir()
    base = ["-run_lengths", "5", "-TR", "1", "-polort", "0", "-num_stimts", "0"]
    output = ["-x1D", "X.xmat.1D", "-x1D_uncensored", "full.xmat.1D"]
    failures = [
        (base + ["-stim_times", "1", "events", "GAM"] + output, 2),
        (["-run_lengths", "0", "-TR", "1", "-polort", "0", "-num_stimts", "0"] + output, 2),
        (["-run_lengths", "5", "-TR", "NaN", "-polort", "0", "-num_stimts", "0"] + output, 2),
        (base + ["-polort", "1"] + output, 2),
        (base + ["-ortvec", "missing.1D", "motion"] + output, 1),
        (base + ["-ortvec", "1D: 1 2", "motion"] + output, 2),
        (base + ["-ortvec", "valid.1D[1]", "motion"] + output, 2),
        (base + ["-ortvec", "valid.1D", "motion", "-TR", "2"] + output, 2),
        (base + ["-ortvec", "valid.1D", "bad;label"] + output, 2),
        (base + ["-ortvec", "valid.1D", "L" * 60] + output, 2),
        (base + ["-ortvec", "valid.1D", "a b", "-ortvec", "valid.1D", "a_b"] + output, 2),
    ]
    (work / "valid.1D").write_text("0\n1\n2\n3\n4\n")
    for name, content in {
        "short": "0\n1\n", "long": "1\n" * 6, "ragged": "1 2\n3\n",
        "nan": "0\n1\nNaN\n3\n4\n", "infinite": "0\n1\n1e999\n3\n4\n",
        "text": "0\n1\noops\n3\n4\n", "empty": "# nothing\n",
    }.items():
        (work / f"{name}.1D").write_text(content)
        failures.append((base + ["-ortvec", f"{name}.1D", "motion"] + output, 2))
    for name, content in {"censor_empty": "0\n" * 5, "censor_fraction": "1\n1\n0.5\n1\n1\n",
                          "censor_columns": "1 1\n" * 5, "censor_short": "1\n"}.items():
        (work / f"{name}.1D").write_text(content)
        failures.append((base + ["-censor", f"{name}.1D"] + output, 2))
    # M4 accepts signal stim_file; missing labels still fail explicitly.
    failures.append((base[:-1] + ["1", "-stim_file", "1", "valid.1D"] + output, 2))
    for args, code in failures:
        invoke(exe, args, work, code)
        assert not (work / "X.xmat.1D").exists()
        assert not (work / "full.xmat.1D").exists()
        assert not list(work.glob(".1dXmatMaker-*.tmp"))

    # Preflight both paths before publishing either, including filesystem aliases.
    (work / "existing").write_text("keep me")
    invoke(exe, base + ["-x1D", "X.xmat.1D", "-x1D_uncensored", "existing"], work, 2)
    assert (work / "existing").read_text() == "keep me"
    assert not (work / "X.xmat.1D").exists()
    invoke(exe, base + ["-x1D", "X.xmat.1D", "-x1D_uncensored", "./X.xmat.1D"], work, 2)
    (work / "dangling").symlink_to("missing-target")
    invoke(exe, base + ["-x1D", "dangling"], work, 2)
    assert (work / "dangling").is_symlink()
    (work / "alias").symlink_to(work, target_is_directory=True)
    invoke(exe, base + ["-x1D", "X.xmat.1D", "-x1D_uncensored", "alias/X.xmat.1D"], work, 2)
    assert not (work / "X.xmat.1D").exists()

    # Default is one output only. Rerunning must preserve its exact bytes.
    invoke(exe, base + ["-x1D", "X.xmat.1D"], work)
    before = digest(work / "X.xmat.1D")
    invoke(exe, base + ["-x1D", "X.xmat.1D"], work, 2)
    assert digest(work / "X.xmat.1D") == before
    assert not (work / "full.xmat.1D").exists()
    # Independent small/wide cases that an OLS fitter may refuse to generate.
    invoke(exe, ["-run_lengths", "1", "-TR", "1", "-polort", "0", "-num_stimts", "0", "-x1D", "one"], work)
    assert read_matrix(work / "one")[1] == [[1.0]]
    invoke(exe, ["-run_lengths", "2", "-TR", "1", "-polort", "3", "-num_stimts", "0", "-x1D", "wide"], work)
    assert len(read_matrix(work / "wide")[1][0]) == 4
    return len(failures) + 8


def baselines_main(argv=None):
    """Capture only on request; otherwise check data, metadata, CLI, and readers."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maker", type=Path)
    parser.add_argument("--reader", type=Path)
    parser.add_argument("--capture-reference", type=Path)
    parser.add_argument("--reference-dir", type=Path, default=BASELINE_REFERENCE)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args(argv)
    cases = json.loads(BASELINE_CASES.read_text())
    if args.capture_reference:
        reference = args.capture_reference.resolve(strict=True)
        pinned = json.loads((FIXTURES / "afni/compatibility/reference/manifest.json").read_text())
        assert digest(reference) == pinned["executable_sha256"], "use the pinned M0 reference binary"
        args.reference_dir.mkdir(exist_ok=False)
        records = [run_case(reference, c, args.reference_dir / c["id"])
                   for c in cases if "m0_case" not in c]
        (args.reference_dir / "manifest.json").write_text(json.dumps({
            "executable_sha256": digest(reference), "cases_sha256": digest(BASELINE_CASES),
            "cases": records}, indent=2) + "\n")
        print(f"Captured {len(records)} additional M1 reference cases")
        return
    if not args.maker:
        parser.error("--maker required for normal validation")
    maker = args.maker.resolve(strict=True)
    reader = args.reader.resolve(strict=True) if args.reader else None
    manifest = json.loads((args.reference_dir / "manifest.json").read_text())
    assert manifest["cases_sha256"] == digest(BASELINE_CASES)
    for entry in manifest["cases"]:
        for name, sha in entry["files"].items():
            assert digest(args.reference_dir / entry["id"] / name) == sha
    reports = []
    with tempfile.TemporaryDirectory(prefix="xmat_m1_check_") as temp:
        work = Path(temp)
        for case in cases:
            directory = work / case["id"]
            directory.mkdir()
            for name, content in case["inputs"].items():
                (directory / name).write_text(content)
            output = ["-x1D", "X.xmat.1D", "-x1D_uncensored", "X.full.xmat.1D"]
            invoke(maker, case["maker_args"] + output, directory)
            expected = FIXTURES / "afni/compatibility/reference" / case["m0_case"] if "m0_case" in case else args.reference_dir / case["id"]
            errors = [compare_baselines(directory / name, expected / name)
                      for name in ("X.xmat.1D", "X.full.xmat.1D")]
            attrs, retained = read_matrix(directory / "X.xmat.1D")
            _, full = read_matrix(directory / "X.full.xmat.1D")
            assert retained == [full[i] for i in integers(attrs["GoodList"])]
            if reader:
                env = {k: v for k, v in os.environ.items() if not k.startswith("AFNI_")}
                env.update(ENV)
                command = [str(reader), "-infile", str(directory / "X.xmat.1D"),
                           "-uncensored", str(directory / "X.full.xmat.1D"), "-describe"]
                result = subprocess.run(command, cwd=work, env=env, capture_output=True, text=True, timeout=30)
                assert result.returncode == 0, result.stderr
                # A successful describe alone could hide incorrectly decoded
                # labels. Compare to the oracle through AFNI's own parser too.
                command = [str(reader), "-infile", str(directory / "X.xmat.1D"),
                           "-compare", str(expected / "X.xmat.1D"), "-atol", "2e-6", "-rtol", "2e-5"]
                result = subprocess.run(command, cwd=work, env=env, capture_output=True, text=True, timeout=30)
                assert result.returncode == 0, (case["id"], result.stdout, result.stderr)
            reports.append({"case": case["id"], "max_abs_difference": max(errors), "afni_reader": bool(reader)})
        boundaries = negative_and_boundary_checks(maker, work / "boundary")
    report = {"maker": str(maker), "maker_sha256": digest(maker),
              "reader": str(reader) if reader else None, "cases": reports,
              "boundary_checks": boundaries, "result": "pass"}
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: {len(reports)} matrix pairs, {boundaries} CLI/boundary checks; AFNI reader={bool(reader)}")


# Event matrices, identities and process helpers (M2)
# ------------------------------------------------------------------------


EVENT_CASES = FIXTURES / 'afni/events/cases.json'
EVENT_REFERENCE = FIXTURES / 'afni/events/reference'


def compare_events(actual, expected):
    """Exact semantic metadata; predeclared float parity and strict TENT gate."""
    aa, ar = read_matrix(actual)
    ea, er = read_matrix(expected)
    assert len(ar) == len(er) and len(ar[0]) == len(er[0])
    ignored = {'CommandLine', 'XmatMakerVersion', 'XmatMakerReference'}
    assert aa.keys()-ignored == ea.keys()-ignored, (actual, aa.keys(), ea.keys())
    for key in ea.keys()-ignored:
        if key in {'GoodList', 'RunStart', 'ColumnGroups', 'StimBots', 'StimTops'}:
            assert integers(aa[key]) == integers(ea[key]), (actual,key)
        elif key == 'RowTR':
            # AFNI prints RowTR using %g; compare to the explicitly supplied
            # acquisition value, not the internal float timing approximation.
            assert abs(float(aa[key])-float(ea[key])) < 1e-14
        elif key == 'ni_type':
            assert aa[key].removeprefix('1*') == ea[key].removeprefix('1*')
        else:
            assert aa[key] == ea[key], (actual,key,aa[key],ea[key])
    labels = [x.strip() for x in ea['ColumnLabels'].split(';')]
    tent_columns = set()
    for key, value in ea.items():
        if key.startswith('BasisFormula_') and value.startswith('TENT'):
            lo, hi = map(int,ea['BasisColumns_'+key.split('_')[1]].split(':'))
            tent_columns.update(range(lo,hi+1))
    worst = 0.0
    for row,(av,ev) in enumerate(zip(ar,er)):
        for col,(a,e) in enumerate(zip(av,ev)):
            strict = col in tent_columns or (labels[col].startswith('Run#') and labels[col].endswith(('Pol#0','Pol#1')))
            tol = 1e-7 if strict else 2e-6+2e-5*abs(e)
            assert abs(a-e) <= tol, (actual,row,col,a,e,tol)
            worst = max(worst,abs(a-e))
    return worst


def expected_relations(case, attrs):
    """Derive event IDs/relations from source tokens, not Maker's manifest.

    Float sampling is intentionally absent: these are original seconds and
    identities. Matrix column ranges come from the independently captured AFNI
    header. Global tokens may interleave runs; IDs and columns have different
    specified orders, which this explicitly checks.
    """
    args=case['maker_args']
    lengths=list(map(int,args[args.index('-run_lengths')+1].split(',')))
    tr=float(args[args.index('-TR')+1]); starts=[sum(lengths[:r]) for r in range(len(lengths))]
    global_times='-global_times' in args
    sources=[]
    for i,option in enumerate(args):
        if option not in ('-stim_times','-stim_times_IM'):continue
        stimulus=int(args[i+1])-1; path=args[i+2]; formula=args[i+3]
        nf=2 if formula=='SPMG2' else 1
        duration=0.0
        if formula.startswith('TENT'):
            nf=int(formula.split(',')[-1].rstrip(')'))-(2 if formula.startswith('TENTzero') else 0)
        elif formula.startswith('BLOCK'):
            duration=float(formula.split('(')[1].split(',')[0].rstrip(')'))
        events=[]; local=0
        for line in case['inputs'][path].splitlines():
            content=line.split('#')[0].strip()
            if not content:continue
            for token in content.split():
                if token=='*':continue
                onset=float(token)
                run=max(r for r,start in enumerate(starts) if start*tr<=onset) if global_times else local
                events.append((run,onset-starts[run]*tr if global_times else onset))
            local+=1
        first=int(attrs[f'BasisColumns_{stimulus+1:06d}'].split(':')[0])
        sources.append((stimulus,option,nf,duration,first,events))
    rows=[]; event_id=0
    for stimulus,option,nf,duration,first,events in sorted(sources):
        for original in sorted(range(len(events)),key=lambda i:(events[i][0],i)):
            run,onset=events[original]
            for basis in range(nf):
                column=first+(original*nf if option=='-stim_times_IM' else 0)+basis
                rows.append({'event':event_id,'stimulus':stimulus,'condition':stimulus,'run':run,
                             'onset':onset,'duration':duration,'basis':basis,'column':column,'item':'-'})
            event_id+=1
    return rows, lengths


def check_manifest(case, directory):
    """Check all identities, and exact zero support outside each event's run."""
    attrs,full=read_matrix(directory/'X.full.xmat.1D')
    expected,lengths=expected_relations(case,attrs)
    text=(directory/'events.tsv').read_text()
    assert text.splitlines()[0]=='# 1dXmatMaker events 1'
    rows=list(csv.DictReader(io.StringIO('\n'.join(text.splitlines()[1:])),delimiter='\t'))
    assert len(rows)==len(expected)
    column_runs={}
    for actual,wanted in zip(rows,expected):
        for key,value in wanted.items():
            if key in ('onset','duration'):assert abs(float(actual[key])-value)<1e-12,(key,actual,wanted)
            elif key=='item':assert actual[key]==value
            else:assert int(actual[key])==value,(key,actual,wanted)
        column_runs.setdefault(wanted['column'],set()).add(wanted['run'])
    for column,runs in column_runs.items():
        start=0
        for run,length in enumerate(lengths):
            if run not in runs: assert all(row[column]==0.0 for row in full[start:start+length])
            start+=length
    return len(rows)


def event_boundaries(exe,work):
    """Reject unsupported/ambiguous event inputs without publishing any sidecar."""
    work.mkdir()
    (work/'times').write_text('2 8\n')
    base=['-run_lengths','30','-TR','1','-polort','0','-num_stimts','1']
    source=['-stim_times_IM','1','times','GAM','-stim_label','1','trial']
    outputs=['-x1D','X','-x1D_uncensored','F','-event_manifest','events']
    checks=0
    def failure(args):
        nonlocal checks
        invoke(exe,args+outputs,work,2)
        assert all(not (work/p).exists() for p in ('X','F','events'))
        assert not list(work.glob('.1dXmatMaker-*.tmp'))
        checks+=1
    failure(base+source) # Never guess the timing mode.
    failure(base+['-local_times','-global_times']+source)
    failure(base+['-local_times']+source+['-stim_base','1'])
    failure(base+['-local_times']+source+['-stim_file','1','times'])
    failure(base+['-local_times']+source+['-global_times'])
    for formula in ['GAM(0,1)','GAM(8.6,.547,-2)','BLOCK(0)','BLOCK(2,-1)','SPMG9','SPMG1(-1)',
                    'TENT(0,3,2.5)','TENTzero(0,3,2)','GAM(1,NaN)','EXPR(0,2)2']:
        altered=source.copy();altered[3]=formula
        failure(base+['-local_times']+altered)
    for text in ['*\n','-1\n','NaN\n','30\n','2*3\n','2:3\n','2 *\n','2\n3\n','2\n*\n']:
        (work/'times').write_text(text)
        failure(base+['-local_times']+source)
    (work/'times').write_text('2 8\n')
    # The manifest joins the same preflight/staging transaction as BOTH matrices.
    (work/'occupied').write_text('keep')
    invoke(exe,base+['-local_times']+source+['-x1D','X','-x1D_uncensored','F','-event_manifest','occupied'],work,2)
    assert (work/'occupied').read_text()=='keep' and not (work/'X').exists() and not (work/'F').exists();checks+=1
    invoke(exe,base+['-local_times']+source+['-x1D','X','-event_manifest','./X'],work,2)
    assert not (work/'X').exists();checks+=1
    # Optional manifest works without a full matrix, and no manifest is default.
    invoke(exe,base+['-local_times']+source+['-x1D','X','-event_manifest','events'],work)
    assert (work/'events').exists() and not (work/'F').exists();checks+=1
    invoke(exe,base+['-local_times']+source+['-x1D','default'],work)
    assert not (work/'default.events.tsv').exists();checks+=1
    return checks


def afni_run(command, directory, expected_error=None):
    """Isolate AFNI config and make the exact executable/arguments reviewable."""
    env={k:v for k,v in os.environ.items() if not k.startswith('AFNI_')}
    env.update(ENV)
    env.update({'AFNI_NIFTI_TYPE_WARN':'NO','AFNI_NIFTI_VIEW':'orig',
                'OPENBLAS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'})
    result=subprocess.run([str(x) for x in command],cwd=directory,env=env,capture_output=True,text=True,timeout=60)
    if expected_error:
        assert result.returncode!=0 and expected_error in result.stderr,(command,result.stdout,result.stderr)
    else:
        assert result.returncode==0,(command,result.stdout,result.stderr)
    return result


def events_main(argv=None):
    """Audit captures, then check matrices/identities/errors and requested consumers."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maker',type=Path)
    parser.add_argument('--afni-bin',type=Path)
    parser.add_argument('--report',type=Path)
    parser.add_argument('--capture-reference',type=Path)
    parser.add_argument('--reference-dir',type=Path,default=EVENT_REFERENCE)
    args=parser.parse_args(argv);cases=json.loads(EVENT_CASES.read_text())
    if args.capture_reference:
        exe=args.capture_reference.resolve(strict=True)
        pinned=json.loads((FIXTURES / 'afni/compatibility/reference/manifest.json').read_text())
        assert digest(exe)==pinned['executable_sha256']
        args.reference_dir.mkdir(exist_ok=False)
        records=[run_case(exe,c,args.reference_dir/c['id']) for c in cases if 'm0_case' not in c]
        (args.reference_dir/'manifest.json').write_text(json.dumps({'executable_sha256':digest(exe),'cases_sha256':digest(EVENT_CASES),'cases':records},indent=2)+'\n')
        return
    if not args.maker:parser.error('--maker is required')
    maker=args.maker.resolve(strict=True)
    afni=args.afni_bin.resolve(strict=True) if args.afni_bin else None
    manifest=json.loads((args.reference_dir/'manifest.json').read_text())
    assert manifest['cases_sha256']==digest(EVENT_CASES)
    for record in manifest['cases']:
        for name,sha in record['files'].items():assert digest(args.reference_dir/record['id']/name)==sha
    records=[];consumer=[]
    with tempfile.TemporaryDirectory(prefix='xmat_m2_') as tmp:
        work=Path(tmp)
        for case in cases:
            directory=work/case['id'];directory.mkdir()
            for name,text in case['inputs'].items():(directory/name).write_text(text)
            invoke(maker,case['maker_args']+['-x1D','X.xmat.1D','-x1D_uncensored','X.full.xmat.1D','-event_manifest','events.tsv'],directory)
            expected=FIXTURES / 'afni/compatibility/reference'/case['m0_case'] if 'm0_case' in case else args.reference_dir/case['id']
            error=max(compare_events(directory/name,expected/name) for name in ('X.xmat.1D','X.full.xmat.1D'))
            attrs,retained=read_matrix(directory/'X.xmat.1D');_,full=read_matrix(directory/'X.full.xmat.1D')
            assert retained==[full[i] for i in integers(attrs['GoodList'])]
            relations=check_manifest(case,directory)
            reader_status='not_run'
            if afni:
                for path in (directory,expected):
                    afni_run([afni/'1dXmatTool','-infile',path/'X.xmat.1D','-uncensored',path/'X.full.xmat.1D','-describe'],directory)
                afni_run([afni/'1dXmatTool','-infile',directory/'X.xmat.1D','-compare',expected/'X.xmat.1D','-atol','2e-6','-rtol','2e-5'],directory)
                reader_status='pass'
                if case['id'] in ('event_censored','mixed_sources','im_spmg_pair'):
                    from .consumers import compare_fits
                    consumer.append(compare_fits(afni,directory,expected,case['id'],afni_run))
            records.append({'case':case['id'],'max_abs_difference':error,'manifest_relations':relations,'afni_reader':reader_status})
        boundaries=event_boundaries(maker,work/'events_boundary')
        old_boundaries=negative_and_boundary_checks(maker,work/'m1_boundary')
    report={'maker':str(maker),'maker_sha256':digest(maker),'afni_bin':str(afni) if afni else None,
            'cases':records,'event_boundary_checks':boundaries,'m1_boundary_checks':old_boundaries,
            'consumer_fits':consumer,'result':'pass'}
    if args.report:args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS: {len(records)} event matrix pairs and manifests; {boundaries} event + {old_boundaries} M1 boundaries; {len(consumer)} consumer cases')

