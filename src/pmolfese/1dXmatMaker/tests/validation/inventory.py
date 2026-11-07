# Developed at the NIH
#
# Part of AFNI. See the repository LICENSE.txt for licensing terms and
# disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
# Third-party and adapted portions retain their existing attribution and license
# terms; see THIRD_PARTY.md.

"""Pinned AFNI basis/options/expression source inventory.

Use tests/run.py for the public commands. Fixtures are grouped by reference
implementation and behavior; checks never rewrite expected numerical results.
"""

import argparse
import hashlib
import json
import re
from . import ROOT
from pathlib import Path



# AFNI source inventory maintenance
# ------------------------------------------------------------------------

REPO = ROOT.parents[2]

# Each design-affecting option has a fixture or a concrete unresolved item.
# Fitting/output options are classified explicitly below instead of silently
# being treated as supported by a matrix-only producer.
DESIGN = {
    "-nodata": "baseline_unequal", "-concat": "baseline_unequal",
    "-polort": "baseline_auto", "-leg": "baseline_quadratic",
    "-legendre": "baseline_quadratic", "-nolegendre": "baseline_power",
    "-dmb": "nuisance", "-dmbase": "nuisance", "-nodmbase": "nuisance_no_demean",
    "-censor": "event_censored", "-CENSOR": "censortr", "-CENSORTR": "censortr",
    "-censorTR": "censortr", "-nfirst": "first_last", "-nlast": "first_last",
    "-num_stimts": "gam", "-stim_file": "nuisance", "-stim_base": "nuisance",
    "-stim_label": "im", "-ortvec": "ortvec", "-inati": "U-ORT-ALIAS",
    "-stim_times": "gam", "-stim_times_IM": "im_multibasis",
    "-stim_times_AM1": "am1", "-stim_times_AM2": "am2",
    "-stim_times_AMx": "amx", "-stim_times_FSL": "fsl", "-stim_times_FS1": "fs1",
    "-local_times": "local_two_runs", "-global_times": "global_two_runs",
    "-guess_times": "U-GUESS", "-stim_times_subtract": "offset_two_runs",
    "-stim_times_millisec": "millisecond_two_runs", "-basis_normall": "basis_norm",
    "-num_glt": "symbolic_glt", "-gltsym": "symbolic_glt", "-glt_label": "symbolic_glt",
    "-glt": "U-GLT-NUMERIC", "-stim_minlag": "U-LAG", "-stim_maxlag": "U-LAG",
    "-stim_nptr": "U-LAG", "-noblock": "U-NOBLOCK",
    "-x1D": "gam", "-x1D_unc": "event_censored", "-x1D_uncensored": "event_censored",
    "-x1D_regcen": "U-REGCEN", "-allzero_OK": "U-ZERO",
    "-goforit": "U-RANK", "-GOFORIT": "U-RANK",
}
OUTSIDE = set('''-OK -notrans -nosat -trans -sat -nocond -sing -singvals -xjpeg
    -nox1D -x1D_stop -input -STATmask -FDRmask -mask -automask -input1D -force_TR
    -TR_1D -noFDR -dofsub -svd -nosvd -xsave -noxsave -xrestore -float -short
    -datum -progress -rmsmin -fdisp -TR_times -TR_irc -slice_base -iresp -tshift
    -sresp -fout -rout -tout -vout -xout -nobout -bout -nocout -full_first
    -nofull_first -nofullf_atall -bucket -prefix -nobucket -cbucket -cprefix
    -fitts -errts -jobs -virtvec -D'''.split())


def generate():
    """Return a classified inventory from local sources and declared cases.

    Unknown option branches or basis families without fixtures fail loudly.
    Special-function presence is inventoried separately from its semantics;
    each function retains an unresolved item until its own tests exist.
    """
    source = (REPO / "src/3dDeconvolve.c").read_text()
    cases = json.loads((ROOT / "fixtures/afni/compatibility/cases.json").read_text())
    # Restrict basis discovery to the actual parser, not help examples.
    begin = source.index("basis_expansion * basis_parser( char *sym )\n{")
    end = source.index("/*--- NO MORE BASIS FUNCTION CHOICES ---*/", begin)
    dispatch = []
    for match in re.finditer(r'(strcmp|strncmp)\(scp,"([^"]+)"', source[begin:end]):
        name = match.group(2)
        if name not in [v["name"] for v in dispatch]:
            dispatch.append({"name": name, "match": match.group(1),
                             "line": source[:begin + match.start()].count("\n") + 1})
    basis = []
    for entry in dispatch:
        name = entry["name"]
        fixtures = []
        for case in cases:
            for i, value in enumerate(case["args"]):
                if value.startswith("-stim_times") and value not in {
                    "-stim_times_millisec", "-stim_times_subtract"}:
                    model = case["args"][i + 3].split("(")[0]
                    if model == name or (entry["match"] == "strncmp" and model.startswith(name)):
                        fixtures.append(case["id"])
        assert fixtures, name
        basis.append({**entry, "implementation": "planned", "fixtures": fixtures,
                      "remaining": "U-PARAMETERS; M3/M4 interactions"})
    options = {}
    pattern = r'(?:strcmp|strncmp|strcasecmp|strstr)\(\s*argv\[nopt\]\s*,\s*"([^"]+)"'
    for match in re.finditer(pattern, source):
        if match.group(1).startswith("-"):
            options[match.group(1)] = source[:match.start()].count("\n") + 1
    for name in set(DESIGN) | OUTSIDE:
        options.setdefault(name, None)
    option_rows = []
    for name, line in sorted(options.items()):
        if name in DESIGN:
            evidence = DESIGN[name]
            status = "unresolved" if evidence.startswith("U-") else "reference-captured"
        elif name in OUTSIDE:
            status, evidence = "outside-core", "See COMPATIBILITY.md scope/exclusions"
        elif name in {"-h", "-help", "-quiet", "-verb"}:
            status, evidence = "interface-only", "M1 help and error presentation"
        else:
            raise AssertionError(f"Unclassified option {name}")
        option_rows.append(dict(option=name, source_line=line, status=status, evidence=evidence))
    # The authoritative parser.inc table includes functions absent from help.
    inc = (REPO / "src/parser.inc").read_text()
    table = inc.split("DATA C_FUNCNAME /", 1)[1].split("'DUMMY'", 1)[0]
    functions = re.findall(r"'([A-Z0-9_]+)'", table)
    count = int(re.search(r"NUM_FUNC = (\d+)", inc).group(1))
    assert len(functions) == count
    expressions = [dict(name=name, status="unresolved", evidence=f"U-EXPR-{name}",
                        gate="M3: arity, domain recovery, accuracy; RNG policy for random functions")
                   for name in functions]
    envs = sorted(set(re.findall(r'(?:getenv|yesenv|noenv)\("([^"]+)"', source)))
    return dict(schema=1, source_sha256=hashlib.sha256(source.encode()).hexdigest(),
                basis_dispatch=basis, options=option_rows, expression_functions=expressions,
                source_environment_variables=envs,
                limitations="Direct comparisons only; shared AFNI preprocessing/options and syntax interactions are tracked in U-INPUT and U-LENIENT.")


def main(argv=None):
    """Write the inventory explicitly, or check it without changing any file."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    target = ROOT / "tests/validation/source_inventory.json"
    result = generate()
    content = json.dumps(result, indent=2) + "\n"
    if args.check:
        assert target.read_text() == content, "Source inventory drift; review before updating"
    else:
        target.write_text(content)
    print(f'PASS: {len(result["basis_dispatch"])} basis dispatch tokens, '
          f'{len(result["options"])} options, {len(result["expression_functions"])} expression functions')

