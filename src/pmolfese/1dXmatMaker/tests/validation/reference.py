# Developed at the NIH
#
# Part of AFNI. See the repository LICENSE.txt for licensing terms and
# disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
# Third-party and adapted portions retain their existing attribution and license
# terms; see THIRD_PARTY.md.

"""Immutable AFNI fixture auditing, replay, and explicit reference capture.

Use tests/run.py for the public commands. Fixtures are grouped by reference
implementation and behavior; checks never rewrite expected numerical results.
"""

import argparse
import hashlib
import html
import json
import math
import os
import re
import subprocess
import tempfile
from . import ROOT
from pathlib import Path



# M0 oracle and independent relations
# ------------------------------------------------------------------------

REPO = ROOT.parents[2]
CASES = ROOT / "fixtures/afni/compatibility/cases.json"
GOLDEN = ROOT / "fixtures/afni/compatibility/reference"

# Clear inherited AFNI settings; -skip_afnirc suppresses the user's .afnirc.
# This is a
# reference policy, not a proposal for how Maker should read its environment.
ENV = {
    "AFNI_NOSPLASH": "YES",
    "AFNI_VERSION_CHECK": "NO", "AFNI_3dDeconvolve_NIML": "YES",
    "AFNI_3dDeconvolve_rawAM2": "NO", "AFNI_3dDeconvolve_extend": "NO",
    "AFNI_3dDeconvolve_GOFORIT": "NO", "AFNI_3dDeconvolve_nodup": "NO",
    "AFNI_3dDeconvolve_nodata_extras": "NO", "LC_ALL": "C",
    "OMP_NUM_THREADS": "1",
}


def digest(path):
    """Hash bytes, not timestamps: rebuilt files may retain identical behavior."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_matrix(path):
    """Read the small reference subset, not a general NIML/AFNI input parser.

    Keep header values as strings, including AFNI's compressed integer lists.
    Semantic decoding is confined to the assertions that need those lists.
    """
    text = Path(path).read_text()
    attrs = dict(re.findall(r'^#\s+(\w+)\s*=\s*"([^"\n]*)"', text, re.M))
    attrs = {key: html.unescape(value) for key, value in attrs.items()}
    rows = [[float(v) for v in line.split()] for line in text.splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    assert rows and len({len(row) for row in rows}) == 1, path
    assert int(attrs["ni_dimen"]) == len(rows), path
    assert all(math.isfinite(v) for row in rows for v in row), path
    return attrs, rows


def integers(text):
    """Decode only the list forms present in this corpus (inclusive ranges)."""
    values = []
    for token in text.split(","):
        if "@" in token:
            count, value = map(int, token.split("@"))
            values.extend([value] * count)
        elif ".." in token:
            first, last = map(int, token.split(".."))
            values.extend(range(first, last + 1))
        else:
            values.append(int(token))
    return values


def check_matrix_pair(directory):
    """Censoring must select rows without changing the underlying regressors."""
    attrs, rows = read_matrix(directory / "X.xmat.1D")
    full_attrs, full = read_matrix(directory / "X.full.xmat.1D")
    good = integers(attrs["GoodList"])
    assert len(good) == len(rows)
    assert good == sorted(set(good))
    assert int(attrs["NRowFull"]) == len(full)
    assert rows == [full[i] for i in good]
    assert len(integers(attrs["ColumnGroups"])) == len(rows[0])
    assert len(attrs["ColumnLabels"].split(";")) == len(rows[0])
    for key in attrs.keys() | full_attrs.keys():
        # The reference omits GLT metadata from its full-row companion.
        # Preserve that fact in the fixtures, rather than inventing metadata.
        per_output = {"ni_dimen", "GoodList", "CommandLine", "Nglt", "GltLabels"}
        if key not in per_output and not key.startswith("GltMatrix_"):
            assert attrs.get(key) == full_attrs.get(key), key


def run_case(executable, case, directory):
    """Materialize declared inputs and execute the exact argv list, no shell."""
    directory.mkdir()
    for name, content in case.get("inputs", {}).items():
        (directory / name).write_text(content)
    env = {k: v for k, v in os.environ.items() if not k.startswith("AFNI_")}
    env.update(ENV)
    env.update(case.get("environment", {}))
    command = [str(executable), "-skip_afnirc"] + case["args"] + [
        "-x1D", "X.xmat.1D", "-x1D_uncensored", "X.full.xmat.1D",
        "-x1D_stop"]
    result = subprocess.run(command, cwd=directory, env=env, text=True,
                            capture_output=True, timeout=60)
    (directory / "stdout.txt").write_text(result.stdout)
    (directory / "stderr.txt").write_text(result.stderr)
    if case["outcome"] == "success":
        if result.returncode != 0:
            raise RuntimeError(f'{case["id"]}: {result.stderr[-2500:]}')
        check_matrix_pair(directory)
    else:
        assert result.returncode != 0, f'{case["id"]}: unexpected success'
        assert case["error_contains"] in result.stderr, case["id"]
    return {"id": case["id"], "returncode": result.returncode,
            "argv": command, "environment": {**ENV, **case.get("environment", {})},
            "files": {p.name: digest(p) for p in sorted(directory.iterdir())
                      if p.is_file()}}


def compare_matrix(actual, expected):
    """Replay tolerances are deliberately tighter than the future Rust port's.

    Command provenance may differ across installations. All other header
    attributes and dimensions must match. This checks numerical repeatability,
    not bitwise reproducibility of platform-dependent math-library functions.
    """
    aa, ar = read_matrix(actual)
    ea, er = read_matrix(expected)
    aa.pop("CommandLine", None)
    ea.pop("CommandLine", None)
    assert aa == ea, f"metadata differs: {actual}"
    assert len(ar) == len(er) and len(ar[0]) == len(er[0])
    for row, (av, ev) in enumerate(zip(ar, er)):
        for col, (a, e) in enumerate(zip(av, ev)):
            assert math.isclose(a, e, abs_tol=1e-12, rel_tol=1e-10), (actual, row, col, a, e)


def independent_checks(corpus):
    """Hand-derived checks and relationships catch errors in fixture selection.

    These checks do NOT establish correctness of every AFNI waveform. They
    anchor a few simple mathematical facts without copying AFNI algorithms.
    """
    def rows(name):
        return read_matrix(corpus / name / "X.full.xmat.1D")[1]

    baseline = rows("baseline_unequal")
    for run, (start, stop) in enumerate(((0, 20), (20, 48))):
        for i in range(start, stop):
            expected = [0.0] * 4
            expected[2 * run] = 1.0
            expected[2 * run + 1] = 2 * (i - start) / (stop - start - 1) - 1
            assert all(math.isclose(a, b, abs_tol=1e-7)
                       for a, b in zip(baseline[i], expected))
    # A single impulse at t=2 and TENT(0,4,3) is three triangular hats.
    tent = rows("analytic_tent")
    for i, row in enumerate(tent):
        expected = [max(0.0, 1 - abs((i - 2) - center) / 2)
                    if 0 <= i - 2 <= 4 else 0.0 for center in (0, 2, 4)]
        assert all(math.isclose(a, b, abs_tol=1e-7) for a, b in zip(row, expected))
    # t is seconds since onset; x and z scale that same time to [0,1], [-1,1].
    for i, row in enumerate(rows("expr_coordinates")):
        t = i - 2
        expected = [t, t / 4, t / 2 - 1] if 0 <= t <= 4 else [0, 0, 0]
        assert all(math.isclose(a, b, abs_tol=1e-7) for a, b in zip(row, expected))
    for a, b in (("local_two_runs", "global_two_runs"),
                 ("local_two_runs", "millisecond_two_runs"),
                 ("local_two_runs", "offset_two_runs"),
                 ("block4", "igfun4"), ("block5", "igfun5"),
                 ("spmg", "spmg2"), ("am2", "amx")):
        av, bv = rows(a), rows(b)
        # Millisecond conversion multiplies floats by 0.001 in AFNI. The
        # nominally equivalent seconds can therefore differ by roundoff.
        assert len(av) == len(bv) and len(av[0]) == len(bv[0]), (a, b)
        assert all(math.isclose(x, y, abs_tol=2e-6, rel_tol=2e-6)
                   for r, s in zip(av, bv) for x, y in zip(r, s)), (a, b)
    centered, raw = rows("am2"), rows("am2_raw")
    # Amplitudes 1,3,5 have mean 3. Baseline is column zero, then mean HRF
    # and modulation: raw modulation = centered modulation + 3*mean HRF.
    for c, r in zip(centered, raw):
        assert math.isclose(c[1], r[1], abs_tol=1e-7)
        assert math.isclose(r[2], c[2] + 3 * c[1], abs_tol=2e-6)
    # A late event in run 1 must not leak into the empty second run.
    for row in rows("edge_empty_run")[20:]:
        assert row[-1] == 0.0
    positive, negative = rows("mion"), rows("mionn")
    for a, b in zip(positive, negative):
        assert math.isclose(a[-1], -b[-1], abs_tol=1e-7)


def main(argv=None):
    """Dispatch explicit capture or read-only checks; fail on mismatched evidence."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["capture", "verify", "audit"])
    parser.add_argument("--executable", type=Path)
    parser.add_argument("--output", type=Path, default=GOLDEN)
    args = parser.parse_args(argv)
    cases = json.loads(CASES.read_text())
    assert len({case["id"] for case in cases}) == len(cases)
    if args.mode == "audit":
        manifest = json.loads((args.output / "manifest.json").read_text())
        assert digest(CASES) == manifest["cases_sha256"]
        assert [c["id"] for c in cases] == [c["id"] for c in manifest["cases"]]
        for record in manifest["cases"]:
            for name, sha in record["files"].items():
                assert digest(args.output / record["id"] / name) == sha
        independent_checks(args.output)
        print(f"PASS: {len(cases)} fixture records, checksums, independent checks")
        return
    if args.executable is None:
        parser.error("capture/verify require --executable")
    exe = args.executable.resolve(strict=True)
    if args.mode == "capture":
        # Refuse replacement; regeneration must be an explicit reviewed change.
        args.output.mkdir(parents=True, exist_ok=False)
        records = [run_case(exe, case, args.output / case["id"]) for case in cases]
        independent_checks(args.output)
        manifest = {
            "schema": 1, "reference_date": "2026-11-06",
            "git_commit": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
            "executable": str(exe), "executable_sha256": digest(exe),
            "cases_sha256": digest(CASES), "cases": records,
        }
        (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    else:
        with tempfile.TemporaryDirectory(prefix="xmat_m0_replay_") as temp:
            for case in cases:
                out = Path(temp) / case["id"]
                run_case(exe, case, out)
                if case["outcome"] == "success":
                    for name in ("X.xmat.1D", "X.full.xmat.1D"):
                        compare_matrix(out / name, args.output / case["id"] / name)
            independent_checks(Path(temp))
    print(f"PASS: {args.mode}, {len(cases)} reference cases")

