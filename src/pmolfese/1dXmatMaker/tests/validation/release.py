# Developed at the NIH
#
# Part of AFNI. See the repository LICENSE.txt for licensing terms and
# disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
# Third-party and adapted portions retain their existing attribution and license
# terms; see THIRD_PARTY.md.

"""Matched performance measurements and offline source-package checks.

Use tests/run.py for the public commands. Fixtures are grouped by reference
implementation and behavior; checks never rewrite expected numerical results.
"""

import argparse
import json
import os
import platform
import statistics
import subprocess
import tempfile
import time
import hashlib
import tarfile
import sys
from . import ROOT
from pathlib import Path
from .reference import read_matrix, digest, integers



# Generation and candidate-batch benchmarks
# ------------------------------------------------------------------------


def measure(command,work):
    """Wait once, collecting per-process peak RSS and monotonic wall seconds."""
    env={k:v for k,v in os.environ.items() if not k.startswith('AFNI_')}
    env.update(OMP_NUM_THREADS='1',AFNI_NIFTI_TYPE_WARN='NO')
    with tempfile.TemporaryFile() as log:
        start=time.perf_counter()
        child=subprocess.Popen(list(map(str,command)),cwd=work,env=env,stdout=log,stderr=log)
        _,status,usage=os.wait4(child.pid,0)
        child.returncode=os.waitstatus_to_exitcode(status)
        elapsed=time.perf_counter()-start
        log.seek(0);message=log.read().decode(errors='replace')
    if child.returncode:raise RuntimeError(f'{command}\n{message}')
    rss=usage.ru_maxrss*(1 if platform.system()=='Darwin' else 1024)
    return dict(wall_seconds=elapsed,peak_rss_bytes=rss)


def summarize(samples):
    return dict(samples=samples,median_wall_seconds=statistics.median(x['wall_seconds'] for x in samples),max_peak_rss_bytes=max(x['peak_rss_bytes'] for x in samples))


def benchmark_main(argv=None):
    """Measure matched generation, separate diagnostics, and candidate batch costs."""
    # NumPy is needed for benchmark parity checks, not offline package checks.
    import numpy as np
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--maker',type=Path,required=True);p.add_argument('--afni-bin',type=Path,required=True)
    p.add_argument('--afni-build',type=Path,help='optional CMake build directory for compiler/flag provenance')
    p.add_argument('--report',type=Path,required=True);p.add_argument('--repeats',type=int,default=3)
    a=p.parse_args(argv);exe=a.maker.resolve(strict=True);afni=a.afni_bin.resolve(strict=True)
    if a.repeats<1:p.error('repeats must be positive')
    report=dict(platform=platform.platform(),machine=platform.machine(),threads=1,repeats=a.repeats,maker_sha256=digest(exe),oracle_sha256=digest(afni/'3dDeconvolve'),workloads=[])
    report['rustc']=subprocess.check_output(['rustc','-Vv'],text=True)
    if a.afni_build:
        cache=(a.afni_build/'CMakeCache.txt').read_text().splitlines()
        report['afni_build_settings']=[line for line in cache if line.startswith(('CMAKE_BUILD_TYPE:','CMAKE_C_COMPILER:','CMAKE_C_FLAGS:','CMAKE_C_FLAGS_RELEASE:','CMAKE_OSX_DEPLOYMENT_TARGET:'))]
        compiler=next(line.split('=',1)[1] for line in cache if line.startswith('CMAKE_C_COMPILER:'))
        report['c_compiler']=subprocess.check_output([compiler,'--version'],text=True)
    with tempfile.TemporaryDirectory(prefix='maker_bench_') as directory:
        root=Path(directory)
        for nrun,nevents in [(200,10),(500,50)]:
            n=2*nrun;record=dict(rows=n,trials=2*nevents,full_row_companion=True,measurements={})
            for mode in ['maker','oracle_no_condition','oracle_default','reader_diagnostics']:
                samples=[]
                for repeat in range(a.repeats):
                    work=root/f'{nrun}_{mode}_{repeat}';work.mkdir()
                    (work/'times').write_text((' '.join(str(3+i*(nrun-20)//nevents) for i in range(nevents))+'\n')*2)
                    (work/'runs').write_text(f'0 {nrun}\n')
                    common=['-nodata',str(n),'1','-concat','runs','-polort','1','-num_stimts','1','-local_times','-stim_times_IM','1','times','GAM','-stim_label','1','trial','-x1D','X.xmat.1D','-x1D_uncensored','X.full.xmat.1D']
                    command=[exe]+common if mode in ['maker','reader_diagnostics'] else [afni/'3dDeconvolve']+common+['-x1D_stop']+(['-nocond'] if mode=='oracle_no_condition' else [])
                    if mode=='reader_diagnostics':
                        # Generation is outside this timed interval: the reader's
                        # full/per-run SVD diagnostics are an optional separate cost.
                        measure(command,work)
                        command=[afni/'1dXmatTool','-infile','X.xmat.1D','-describe']
                    samples.append(measure(command,work))
                    current=read_matrix(work/'X.xmat.1D');full=read_matrix(work/'X.full.xmat.1D')
                    if mode=='maker' and repeat==0:expected=current;expected_full=full
                    for actual,reference in [(current,expected),(full,expected_full)]:
                        np.testing.assert_allclose(actual[1],reference[1],atol=1e-6,rtol=1e-6)
                        for field in ['ColumnLabels','ColumnGroups','RunStart','GoodList','StimBots','StimTops']:
                            left,right=actual[0].get(field),reference[0].get(field)
                            if field!='ColumnLabels' and left is not None and right is not None:left,right=integers(left),integers(right)
                            assert left==right,field
                record['measurements'][mode]=summarize(samples)
            record['numeric_parity']='pass';report['workloads'].append(record)
        report['candidate_batches']=[]
        for count in [2,20]:
            samples=[]
            for repeat in range(a.repeats):
                work=root/f'candidate_{count}_{repeat}';work.mkdir()
                (work/'events').write_text('trial\tcondition\trun\tonset\tduration\titem\n'+''.join(f'{i}\t{i%4}\t{i//25}\t{3+18*(i%25)}\t1\t-\n' for i in range(50)))
                command=[exe,'-candidate_library','glmsingle','-candidate_ids',','.join(map(str,range(count))),'-events','events','-run_lengths','500,500','-TR','1','-polort','1','-candidate_prefix','hrf']
                samples.append(measure(command,work))
            outputs=list(work.glob('hrf*'))
            report['candidate_batches'].append(dict(candidates=count,rows=1000,trials=50,files=len(outputs),bytes=sum(p.stat().st_size for p in outputs),measurement=summarize(samples)))
    report['notes']=['Maker performs finite/input checks but no condition-number SVD. Oracle -nocond is the closest generation-only baseline; defaults are recorded separately.','Reader diagnostics have their own algorithms and are not equivalent to the oracle condition checks.','Candidate batches include trial and ONOFF matrices, retained and full copies, waveforms and identities. No C baseline is claimed for GLMsingle.','One full design is cached; full/retained output share it. HRF curves and publication paths scale with candidate count.','Times include process startup and text/file output, with warm filesystem caches; these are synthetic native measurements, not universal speed claims.']
    a.report.write_text(json.dumps(report,indent=2)+'\n');print('PASS matched numeric benchmark workloads and candidate batches')


# Archive, offline build and reproducibility checks
# ------------------------------------------------------------------------


def run_command(command,work,env):
    """Retain useful failure diagnostics without flooding a successful report."""
    result=subprocess.run(list(map(str,command)),cwd=work,env=env,capture_output=True,text=True)
    if result.returncode:raise RuntimeError(f'{command}\n{result.stdout}\n{result.stderr}')
    return result.stdout+result.stderr


def package_main(argv=None):
    """Verify archives and build/test their real source tree with an empty Cargo cache."""
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--release-dir',type=Path,required=True);p.add_argument('--report',type=Path,required=True);a=p.parse_args(argv)
    release=a.release_dir.resolve(strict=True);hashes=json.loads((release/'SHA256SUMS.json').read_text())
    for name,expected in hashes.items():
        if Path(name).name!=name:raise ValueError('checksum path must be a basename')
        assert hashlib.sha256((release/name).read_bytes()).hexdigest()==expected,name
    source=next(name for name in hashes if name.endswith('-source.tar.gz'))
    native=next(name for name in hashes if name!=source)
    with tempfile.TemporaryDirectory(prefix='maker_package_') as directory:
        root=Path(directory)
        for kind,name in [('source',source),('native',native)]:
            destination=root/kind;destination.mkdir()
            with tarfile.open(release/name) as archive:archive.extractall(destination,filter='data')
        project=next((root/'source').iterdir());distribution=next((root/'native').iterdir())
        env={k:v for k,v in os.environ.items() if not k.startswith(('AFNI_','CARGO_','RUSTFLAGS','CARGO_ENCODED_RUSTFLAGS'))}
        cache=root/'empty-cargo-home';cache.mkdir();env['CARGO_HOME']=str(cache);env['CARGO_TARGET_DIR']=str(root/'build')
        build=run_command(['cargo','build','--release','--bins','--examples','--frozen'],project,env)
        tests=run_command(['cargo','test','--frozen'],project,env)
        # Validate the relocated test entry point and fixtures from the extracted
        # release itself, without relying on this checkout or scientific Python.
        run_command([sys.executable,project/'tests/run.py','--help'],project,env)
        audit=run_command([sys.executable,project/'tests/run.py','reference','audit'],project,env)
        executables=[distribution/'bin/1dXmatMaker',root/'build/release/1dXmatMaker'];sets=[]
        for i,exe in enumerate(executables):
            work=root/f'example{i}';work.mkdir()
            # Identical argv matter for the command-history header. All example
            # inputs are copied under the same relative names in each cwd.
            import shutil
            shutil.copytree(project/'examples',work/'examples')
            run_command([exe,'-candidate_library','glmsingle','-candidate_ids','0,10,19','-events','examples/candidate_events.tsv','-run_lengths','144,144,144','-TR','1','-polort','0','-candidate_prefix','bank'],work,env)
            run_command([exe,'-candidate_library','examples/custom_library.tsv','-hrf_normalization','none','-events','examples/candidate_events.tsv','-run_lengths','144,144,144','-TR','1','-polort','0','-candidate_prefix','custom'],work,env)
            run_command([exe,'-run_lengths','20,28','-TR','1','-polort','1','-num_stimts','0','-x1D','ordinary'],work,env)
            sets.append({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in work.iterdir() if p.is_file()})
        assert sets[0]==sets[1],'packaged/source-built example outputs differ'
        report=dict(result='pass',archives=hashes,offline_build='cargo build --release --bins --examples --frozen with empty CARGO_HOME',packaged_runner='help and fixture audit pass',packaged_probes='basis_probe and hrf_probe built',unit_tests=[line for line in tests.splitlines() if line.startswith('test result:')],example_files=len(sets[0]),example_reproducibility='byte-identical packaged and source-built outputs',packaged_binary_sha256=hashlib.sha256(executables[0].read_bytes()).hexdigest(),rebuilt_binary_sha256=hashlib.sha256(executables[1].read_bytes()).hexdigest())
    a.report.write_text(json.dumps(report,indent=2)+'\n');print('PASS archive checksums, empty-cache frozen build/tests, and matching examples')

