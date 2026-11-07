#!/usr/bin/env python3
# Developed at the NIH
#
# Part of AFNI. See the repository LICENSE.txt for licensing terms and
# disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
# Third-party and adapted portions retain their existing attribution and license
# terms; see THIRD_PARTY.md.

"""One entry point for scientific checks, reference maintenance and release QA.

Run `python3 tests/run.py suite --help` for the complete scientific gate, or
`python3 tests/run.py COMMAND --help` for a focused check. Normal suite runs
never capture expected values. Rust unit/integration tests remain `cargo test`.
"""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import subprocess
import sys

from validation import ROOT

# Keep descriptions and routing in one table. Modules are imported lazily so
# --help, fixture audits, and archive checks do not require NumPy or SciPy.
COMMANDS = {
    'baselines': ('designs', 'baselines_main', 'baseline/nuisance matrices and CLI boundaries'),
    'events': ('designs', 'events_main', 'event matrices, identities and basic consumer fits'),
    'basis': ('basis', 'basis_main', 'frozen basis/EXPR oracle and grammar checks'),
    'math': ('basis', 'math_main', 'independent SciPy numerical anchors'),
    'semantics': ('semantics', 'semantics_main', 'timing, modulation, censoring and contrasts'),
    'contracts': ('semantics', 'contracts_main', 'rejection, publication and environment checks'),
    'consumers': ('consumers', 'consumers_main', 'REML/Picard/LSS, run splits and censor spikes'),
    'sampled': ('sampled', 'sampled_main', 'pinned GLMsingle and custom sampled HRFs'),
    'workflow': ('consumers', 'workflow_main', 'candidate selection and RSA/MVPA identities'),
    'benchmark': ('release', 'benchmark_main', 'matched generation timings and per-child RSS'),
    'package': ('release', 'package_main', 'archive checksums, offline build and examples'),
    'reference': ('reference', 'main', 'audit/replay AFNI compatibility fixtures or explicitly capture'),
    'inventory': ('inventory', 'main', 'check or explicitly update the AFNI source inventory'),
    'capture-hrfs': ('sampled', 'capture_main', 'explicit upstream capture into a NEW fixture file'),
}


def digest(path):
    """Hash the exact tested artifact or report, without importing scientific code."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def suite(argv):
    """Run each gate in a child process, retaining diagnostics and immutable goldens.

    Separate processes prevent one check's argparse/global state from changing
    another's defaults. We use this same runner for children: no obsolete script
    paths or hidden alternate CLI remain. Generated reports live under target/.
    """
    parser = argparse.ArgumentParser(prog='tests/run.py suite', description=suite.__doc__)
    parser.add_argument('--maker', type=Path, required=True)
    parser.add_argument('--afni-bin', type=Path, required=True)
    parser.add_argument('--probes', type=Path, required=True,
                        help='directory containing basis_probe and hrf_probe')
    parser.add_argument('--tag', default='release')
    parser.add_argument('--report-dir', type=Path, default=ROOT/'target/validation')
    args = parser.parse_args(argv)
    if not args.tag or any(c not in 'abcdefghijklmnopqrstuvwxyz0123456789_-' for c in args.tag):
        parser.error('tag must contain lowercase letters, digits, underscores or hyphens')
    maker = args.maker.resolve(strict=True)
    afni = args.afni_bin.resolve(strict=True)
    probes = args.probes.resolve(strict=True)
    output = args.report_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    for path in [afni/name for name in ['1dXmatTool','3dPicard','3dREMLfit','3dLSS','3dMVPA','3dRSA','3dinfo']] + [probes/'basis_probe',probes/'hrf_probe']:
        if not path.is_file():
            parser.error(f'missing required test executable: {path}')
    maker_args = ['--maker', maker]
    cases = [
        ('baselines', 'baselines', maker_args+['--reader', afni/'1dXmatTool']),
        ('events', 'events', maker_args+['--afni-bin', afni]),
        ('basis', 'basis', maker_args+['--probe', probes/'basis_probe','--reader', afni/'1dXmatTool']),
        ('math', 'math', ['--probe', probes/'basis_probe']),
        ('semantics', 'semantics', maker_args+['--reader', afni/'1dXmatTool']),
        ('contracts', 'contracts', maker_args),
        ('consumers', 'consumers', maker_args+['--afni-bin', afni]),
        ('sampled', 'sampled', maker_args+['--probe', probes/'hrf_probe','--picard', afni/'3dPicard']),
        ('workflow', 'workflow', maker_args+['--afni-bin', afni]),
    ]
    reports = {}
    # Audit first: test expectations must be intact before running the producer.
    commands = [('audit', 'reference', ['audit'], None)]
    commands += [(name, command, options, output/f'{args.tag}-{name}.json')
                 for name, command, options in cases]
    summary_path = output/f'{args.tag}-suite.json'
    summary = dict(result='running', maker_sha256=digest(maker), reports=reports,
                   test_source_sha256={str(p.relative_to(ROOT)):digest(p)
                                       for p in [Path(__file__).resolve(), *sorted((ROOT/'tests/validation').glob('*.py'))]})
    # Replace a stale success summary immediately, including when a later gate fails.
    summary_path.write_text(json.dumps(summary, indent=2)+'\n')
    for name, command, options, report in commands:
        if report is not None:
            report.unlink(missing_ok=True)
            options = options+['--report', report]
        log = output/f'{args.tag}-{name}.log'
        invocation = [sys.executable, Path(__file__).resolve(), command, *options]
        with log.open('w') as stream:
            result = subprocess.run(list(map(str, invocation)), stdout=stream, stderr=subprocess.STDOUT)
        if result.returncode or (report is not None and not report.is_file()):
            summary.update(result='FAIL', failed=name, log=str(log))
            summary_path.write_text(json.dumps(summary, indent=2)+'\n')
            raise SystemExit(f'{name} failed; see {log}\n{log.read_text()[-6000:]}')
        if report is not None:
            reports[report.name] = digest(report)
        print(f'{name}: PASS', flush=True)
    summary['result'] = 'pass'
    summary_path.write_text(json.dumps(summary, indent=2)+'\n')
    print(f'PASS: fixture audit and {len(cases)} scientific gates; {summary_path}')


def main():
    """Dispatch one explicit command; importing modules never launches a suite."""
    descriptions = ['  suite          all scientific gates, including fixture audit and SciPy']
    descriptions += [f'  {name:14} {entry[2]}' for name, entry in COMMANDS.items()]
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog='Commands:\n'+'\n'.join(descriptions))
    parser.add_argument('command', choices=['suite', *COMMANDS])
    arguments = sys.argv[1:]
    if not arguments or arguments[0] in ('-h','--help'):
        parser.print_help()
        return
    command = parser.parse_args(arguments[:1]).command
    if command == 'suite':
        suite(arguments[1:])
    else:
        module, function, _ = COMMANDS[command]
        getattr(importlib.import_module('validation.'+module), function)(arguments[1:])


if __name__ == '__main__':
    main()
