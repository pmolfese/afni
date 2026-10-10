#!/usr/bin/env python3
"""End-to-end checks of the 3dSRM program (fit mode) on generated data.

Builds a small multi-subject data set (NIfTI volumes with offsets, scales, constant voxels, and an
all-zero volume; per-subject masks; and 1D inputs), runs the real 3dSRM on it, and compares every
output with an independent computation: the library driver (srm_cli, which is the same thd_srm.c
but with none of the table/dataset/preprocessing code) run on arrays preprocessed here in numpy.

It also asserts the logging convention: the program name appears on the banner line only.

  python check_3dSRM.py --prog build/targets_built/3dSRM --workdir ./check_3dsrm

Threads are capped and the programs run at low priority, so it is safe on a working machine.
Requires numpy and nibabel. The run is SKIPPED, not failed, if either is
missing, unless --require-deps is given; CTest gives it, so a missing environment cannot pass.
AFNI_SRM_KERNEL in the environment (VOXEL or GRAM) reaches 3dSRM, so the same checks can be run
on either fitting kernel.
"""

import argparse
import hashlib
import os
import re
import subprocess
import sys
from pathlib import Path

try:
    import nibabel as nib
    import numpy as np
except ImportError as e:
    if "--require-deps" in sys.argv:
        sys.exit(f"missing a required module: {e}")
    print(f"SKIP: {e}")
    sys.exit(0)

REPO = Path(__file__).resolve().parents[2]

ENV = dict(os.environ, OMP_NUM_THREADS="4", OPENBLAS_NUM_THREADS="4", VECLIB_MAXIMUM_THREADS="4",
           MKL_NUM_THREADS="4", NUMEXPR_NUM_THREADS="4")


# ------------------------------------------------------------------ the library driver

def build_cli(workdir, given, openmp=True):
    """Compile srm_cli from the repo. With OpenMP when Homebrew's libomp is
    present: the fit is threaded, and the results do not depend on the threads."""
    if given:
        return Path(given)
    src = REPO / "src" / "pmolfese" / "tests" / "srm_cli.c"
    thd = REPO / "src" / "thd_srm.c"
    omp = Path("/opt/homebrew/opt/libomp")
    use_omp = openmp and (omp / "include" / "omp.h").exists()
    exe = workdir / ("srm_cli_omp" if use_omp else "srm_cli")
    if not exe.exists() or exe.stat().st_mtime < max(src.stat().st_mtime, thd.stat().st_mtime,
                                                     (REPO / "src" / "thd_srm.h").stat().st_mtime):
        cmd = ["cc", "-O2", "-I", str(REPO / "src"), "-o", str(exe), str(src), "-lm"]
        if use_omp:
            cmd[2:2] = ["-DUSE_OMP", "-Xpreprocessor", "-fopenmp", f"-I{omp}/include",
                        f"-L{omp}/lib", "-lomp"]
        subprocess.run(cmd, check=True)
    return exe


def run_ours(cli, workdir, tag, algo, views, k, niter, center, init_w=None,
             nstart=1, tol=0.0, seed=1234567, extra=None, seglen=None, coupling=None):
    """Fit the views (each V_i x T) with the driver; returns S, W, the objective trace, and the
    noise model (prob), as float64 arrays."""
    d = workdir / tag
    d.mkdir(parents=True, exist_ok=True)
    T = views[0].shape[1]
    cmd = [str(cli), "--algo", algo, "--k", str(k), "--niter", str(niter),
           "--center", str(int(center)), "--T", str(T), "--out", str(d / "out"),
           "--nstart", str(nstart), "--tol", repr(tol), "--seed", str(seed)]
    if init_w is not None:
        for i, w in enumerate(init_w):
            w.astype(np.float32).tofile(d / f"init.W{i}.f32")
        cmd += ["--init", str(d / "init")]
    for i, x in enumerate(views):
        f = d / f"x{i}.f32"
        a = np.ascontiguousarray(x, dtype=np.float32)
        stamp = d / f"x{i}.md5"
        digest = hashlib.md5(a.tobytes()).hexdigest()       # a size match alone is not a content match
        if not f.exists() or not stamp.exists() or stamp.read_text() != digest:
            a.tofile(f)
            stamp.write_text(digest)
        cmd += ["--view", str(x.shape[0]), str(f)]
    cmd += list(extra or [])
    if coupling is not None:
        cmd += ["--seglen", ",".join(str(s) for s in seglen), "--coupling", repr(float(coupling))]
    subprocess.run(cmd, check=True)
    out = {}
    out["s"] = np.fromfile(d / "out.S.f32", dtype=np.float32).reshape(k, T).astype(np.float64)
    out["w"] = [np.fromfile(d / f"out.W{i}.f32", dtype=np.float32).reshape(-1, k).astype(np.float64)
                for i in range(len(views))]
    out["obj"] = np.fromfile(d / "out.obj.f64", dtype=np.float64)
    if coupling is not None:
        out["wr"] = [[np.fromfile(d / f"out.Wr{i}_{r}.f32", dtype=np.float32).reshape(-1, k).astype(np.float64)
                      for r in range(len(seglen))] for i in range(len(views))]
    out["starts"] = np.fromfile(d / "out.starts.f64", dtype=np.float64)
    if algo == "prob":
        out["rho2"] = np.fromfile(d / "out.rho2.f64", dtype=np.float64)
        out["sigma"] = np.fromfile(d / "out.sigma.f32", dtype=np.float32).reshape(k, k).astype(np.float64)
    niter_done, converged, _, best = open(d / "out.info.txt").read().split()
    out["niter"], out["converged"], out["best"] = int(niter_done), bool(int(converged)), int(best)
    return out


def run_load_update(cli, model_file, out, new_files, new_rows, niter):
    cmd = [str(cli), "--load-model", str(model_file), "--out", str(out), "--niter", str(niter), "--tol", "0"]
    for f, r in zip(new_files, new_rows):
        cmd += ["--update-view", str(r), str(f)]
    subprocess.run(cmd, check=True)


def parse_srm_file(path):
    """A reader written from the format described in thd_srm.h, not from the C code."""
    hdr, subj, items, blocks = {}, [], [], {}
    lines = Path(path).read_text().splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i].strip()
        i += 1
        if not ln or ln.startswith("#"):
            continue
        key, *rest = ln.split()
        if key == "data":
            name, r, c = rest[0], int(rest[1]), int(rest[2])
            vals = []
            while len(vals) < r * c:
                vals += lines[i].split()
                i += 1
            blocks[name] = np.array(vals[:r * c], dtype=np.float64).reshape(r, c)
        elif key == "subject":
            subj.append((rest[0], int(rest[1])))
        elif key == "item":
            items.append(rest[0])
        elif key == "end":
            break
        else:
            hdr[key] = rest
    return hdr, subj, items, blocks


def relfro(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def procrustes_residual(a, b):
    """|R a - b| / |b| for the best orthogonal R."""
    u, _, vt = np.linalg.svd(b @ a.T)
    return relfro(u @ vt @ a, b)

failures = []


def check(label, value, limit, fmt="{:.2e}"):
    ok = value <= limit
    print(f"    {label:70s} {fmt.format(value):>10s}  (limit {limit:g})  {'ok' if ok else 'FAIL'}")
    if not ok:
        failures.append(label)


def run(cmd, expect_ok=True):
    p = subprocess.run(["nice", "-n", "15"] + [str(x) for x in cmd], capture_output=True, text=True, env=ENV)
    if expect_ok and p.returncode != 0:
        print(p.stdout, p.stderr)
        raise SystemExit(f"command failed: {' '.join(str(x) for x in cmd)}")
    return p


# --------------------------------------------------------------------- data

NX, NY, NZ = 10, 10, 6
RUNS = [60, 50, 70]
K = 3


def make_world(rng, nsub=8):
    T = sum(RUNS)
    s = np.zeros((K, T))
    c0 = 0
    for L in RUNS:                                   # smooth shared signal, restarted in every run
        s[:, c0] = rng.standard_normal(K)
        for t in range(1, L):
            s[:, c0 + t] = 0.7 * s[:, c0 + t - 1] + np.sqrt(1 - 0.49) * rng.standard_normal(K)
        c0 += L
    s *= np.array([[1.6], [1.1], [0.7]])
    subs = []
    for i in range(nsub):
        mask = np.zeros((NX, NY, NZ), dtype=bool)
        if i < nsub // 2:
            mask[1:9, 1:9, 1:5] = True               # 256 voxels
        else:
            mask[2:8, 1:9, 1:5] = True               # 192 voxels: a different mask per subject
        idx = np.flatnonzero(mask.ravel(order="F"))  # AFNI order: x fastest
        v = len(idx)
        q, _ = np.linalg.qr(rng.standard_normal((v, K)))
        sigma = 0.6 + 0.15 * i
        x = q @ s + sigma * rng.standard_normal((v, T))
        x = x * rng.uniform(0.5, 3.0, (v, 1)) + rng.uniform(50, 500, (v, 1))     # scale and offset
        subs.append({"label": f"sub{i + 1:02d}", "group": "ctl" if i % 2 == 0 else "pat", "mask": mask,
                     "idx": idx, "x": x})
    # a voxel that never varies, one that is flat in a single run, and an all-zero volume
    subs[1]["x"][5, :] = 7.0
    if len(RUNS) > 1:
        subs[2]["x"][9, RUNS[0]:RUNS[0] + RUNS[1]] = 3.0
    return s, subs


def write_nifti(path, arr4d):
    nib.save(nib.Nifti1Image(arr4d.astype(np.float32), np.diag([2.0, 2.0, 2.0, 1.0])), str(path))


def write_world(d, subs, zero_volume=None):
    d.mkdir(parents=True, exist_ok=True)
    rows = []
    for sub in subs:
        c0 = 0
        mask_file = d / f"{sub['label']}_mask.nii.gz"
        write_nifti(mask_file, sub["mask"].astype(np.float32)[..., None])
        for r, L in enumerate(RUNS):
            vol = np.zeros((NX, NY, NZ, L), dtype=np.float32)
            flat = vol.reshape(-1, L, order="F")
            flat[sub["idx"], :] = sub["x"][:, c0:c0 + L]
            vol = flat.reshape((NX, NY, NZ, L), order="F")
            if zero_volume and zero_volume[0] == sub["label"] and zero_volume[1] == r:
                vol[..., zero_volume[2]] = 0.0
            f = d / f"{sub['label']}_run{r + 1}.nii.gz"
            write_nifti(f, vol)
            rows.append((sub["label"], sub["group"], f"r{r + 1}", mask_file, f))
            c0 += L
    return rows


def write_table(path, rows, with_mask):
    lines = ["Subj Group Run " + ("Mask " if with_mask else "") + "InputFile"]
    for s, g, r, m, f in rows:
        lines.append(f"{s} {g} {r} " + (f"{m} " if with_mask else "") + f"{f}")
    Path(path).write_text("\n".join(lines) + "\n")


def load_afni(prefix):
    """An output dataset by prefix, whichever form AFNI wrote it in."""
    for suf in ("+orig.HEAD", ".nii.gz", ".nii"):
        f = Path(str(prefix) + suf)
        if f.exists():
            return np.asarray(nib.load(str(f)).dataobj)
    raise FileNotFoundError(prefix)


def preprocess(sub, keep=None, zscore=True):
    """What 3dSRM should do: censor, z-score within runs, drop features that never vary."""
    x = sub["x"].astype(np.float32).astype(np.float64)
    cols = np.arange(x.shape[1]) if keep is None else np.flatnonzero(keep)
    x = x[:, cols]
    lens, c0 = [], 0
    for L in RUNS:
        lens.append(int(((cols >= c0) & (cols < c0 + L)).sum()))
        c0 += L
    z = np.zeros_like(x)
    bad = np.zeros(x.shape[0], dtype=bool)
    c0 = 0
    for L in lens:
        seg = x[:, c0:c0 + L]
        sd = seg.std(axis=1)
        bad |= ~(sd > 1e-7 * np.abs(seg.mean(axis=1)))
        z[:, c0:c0 + L] = (seg - seg.mean(axis=1, keepdims=True)) / np.where(sd > 0, sd, 1.0)[:, None]
        c0 += L
    return z[~bad], ~bad, lens


# --------------------------------------------------------------------- checks

def reference_fit(cli, workdir, tag, arrays, k, niter, algo, extra=None):
    views = [a.astype(np.float32) for a in arrays]
    return run_ours(cli, workdir, tag, algo, views, k, niter, False, extra=["--canon"] + (extra or []))


def case_basic(prog, cli, work, rng, with_mask_col, zcensor, algo):
    name = f"{'maskcol' if with_mask_col else 'globalmask'}_{'zc' if zcensor else 'nozc'}_{algo}"
    print(f"\n[{name}]", flush=True)
    d = work / name
    s_true, subs = make_world(rng)
    zero = ("sub04", 1, 10) if zcensor else None
    rows = write_world(d, subs, zero)
    write_table(d / "table.txt", rows, with_mask_col)
    # global mask: the union geometry must serve every subject, so use the smaller common mask here
    if not with_mask_col:
        common = np.logical_and.reduce([s["mask"] for s in subs])
        write_nifti(d / "mask.nii.gz", common.astype(np.float32)[..., None])
        for s in subs:                                # restrict every subject to the common mask
            sel = common.ravel(order="F")[s["idx"]]
            s["x"], s["idx"] = s["x"][sel], s["idx"][sel]
        rows = write_world(d, subs, zero)
        write_table(d / "table.txt", rows, False)
    niter = 60
    cmd = [prog, "-prefix", d / "out", "-nfeatures", K, "-algorithm", algo, "-niter", niter, "-tol", 0,
           "-dataTableFile", d / "table.txt", "-save_W", "-progress", "line", "-overwrite"]
    if not with_mask_col:
        cmd += ["-mask", d / "mask.nii.gz"]
    if zcensor:
        cmd += ["-zcensor"]
    p = run(cmd)
    log = p.stdout + p.stderr
    # the logging convention: the program is named once, on the banner
    named = [ln for ln in log.splitlines() if "3dSRM" in ln]
    check("the program name appears on the banner line only", float(len(named) - 1), 0.0, "{:.0f}")
    check("every message line starts with '++', '*+', '**' or ' +' (AFNI style)",
          float(sum(1 for ln in log.splitlines() if ln.strip() and not ln.startswith(("++", "*+", "**", " +", "\r")))),
          0.0, "{:.0f}")

    keep = None
    if zcensor:
        keep = np.ones(sum(RUNS), dtype=bool)
        keep[RUNS[0] + 10] = False
    pre = [preprocess(s, keep) for s in subs]
    T = pre[0][0].shape[1]
    ref = reference_fit(cli, work, name + "_ref", [z for z, _, _ in pre], K, niter, algo)

    S = np.loadtxt(d / "out.S.1D").T                     # file: one row per time point
    check("S.1D vs the library driver, relative difference", relfro(S, ref["s"]), 2e-4)
    fit = (d / "out.fit.txt").read_text().splitlines()
    rows_fit = [ln.split("\t") for ln in fit if not ln.startswith("#") and not ln.startswith("Subj")]
    used = [int(r[3]) for r in rows_fit]
    check("features used per subject (constant ones dropped)", float(sum(abs(u - ok.sum()) for u, (_, ok, _) in zip(used, pre))), 0.0, "{:.0f}")
    if algo == "prob":
        rho = np.array([float(r[5]) for r in rows_fit])
        check("rho2 vs the library driver, max relative", float(np.abs(rho - ref["rho2"]).max() / ref["rho2"].max()), 2e-3)
    # projected data and the maps
    worst_p = worst_w = 0.0
    for i, sub in enumerate(subs):
        z, ok, _ = pre[i]
        proj = load_afni(d / f"out.proj.{sub['label']}").reshape(K, T)
        worst_p = max(worst_p, relfro(proj, ref["w"][i].T @ z))
        wv = load_afni(d / f"out.W.{sub['label']}").reshape(-1, K, order="F")
        worst_w = max(worst_w, relfro(wv[sub["idx"][ok]], ref["w"][i]))
    check("projected series (dataset) vs W'X from the driver, worst subject", worst_p, 5e-3)
    check("W maps (-save_W) vs the driver, worst subject", worst_w, 5e-3)
    # the model file
    hdr, subj, _, blocks = parse_srm_file(d / "out.srm")
    check("model file: S matches S.1D", float(np.abs(blocks["S"] - S).max()), 2e-6)
    check("model file: runs and censoring recorded",
          float(int(hdr["nseg"][0]) != 3 or int(hdr["ntime_in"][0]) != sum(RUNS) or int(hdr["ntime"][0]) != T), 0.0, "{:.0f}")
    return log


def case_1D_and_errors(prog, cli, work, rng):
    print("\n[1D inputs, a mixed table, and refusals]", flush=True)
    d = work / "mixed"
    s_true, subs = make_world(rng, nsub=6)
    d.mkdir(parents=True, exist_ok=True)
    # two "EEG" subjects as 1D files, features = channels; they hear the same latent signal
    eeg = []
    for e in range(2):
        nch = 10 + 4 * e
        q, _ = np.linalg.qr(rng.standard_normal((nch, K)))
        x = (q @ s_true + 0.9 * rng.standard_normal((nch, sum(RUNS)))) * rng.uniform(0.5, 2, (nch, 1)) + 20
        eeg.append({"label": f"eeg{e + 1}", "group": "ctl", "x": x})
    rows = write_world(d, subs)
    lines = ["Subj Group Run Mask InputFile"] + [f"{s} {g} {r} {m} {f}" for s, g, r, m, f in rows]
    for e in eeg:
        c0 = 0
        for r, L in enumerate(RUNS):
            f = d / f"{e['label']}_run{r + 1}.1D"
            np.savetxt(f, e["x"][:, c0:c0 + L].T, fmt="%.8g")        # time x channels
            lines.append(f"{e['label']} {e['group']} r{r + 1} - {f}")
            c0 += L
    (d / "table.txt").write_text("\n".join(lines) + "\n")
    out = d / "out"
    p = run([prog, "-prefix", out, "-nfeatures", K, "-niter", 40, "-tol", 0, "-dataTableFile", d / "table.txt",
             "-progress", "off", "-save_W", "-proj_format", "1D", "-overwrite"])
    names = [s["label"] for s in subs] + [e["label"] for e in eeg]
    pre = [preprocess(s)[0] for s in subs] + [preprocess({"x": e["x"]})[0] for e in eeg]
    ref = reference_fit(cli, work, "mixed_ref", pre, K, 40, "prob")
    S = np.loadtxt(d / "out.S.1D").T
    check("a table mixing datasets and 1D subjects: S vs the driver", relfro(S, ref["s"]), 2e-4)
    pe = np.loadtxt(d / "out.proj.eeg1.1D").T
    check("a 1D subject's projected series (.1D) vs W'X", relfro(pe, ref["w"][len(subs)].T @ pre[len(subs)]), 5e-3)
    wmat = np.loadtxt(d / "out.W.eeg2.1D")
    check("a 1D subject's W (.1D), one row per channel", relfro(wmat, ref["w"][len(subs) + 1]), 5e-3)

    # refusals, each with a message that says what is wrong
    def refuses(label, args, text, table=None):
        t = table or (d / "table.txt")
        pp = run([prog, "-prefix", d / "ref", "-nfeatures", K, "-progress", "off", "-dataTableFile", t, "-overwrite"] + args,
                 expect_ok=False)
        check(f"refuses: {label}", float(not (pp.returncode != 0 and text in (pp.stdout + pp.stderr))), 0.0, "{:.0f}")

    refuses("-cv_adapt (not in this build)", ["-cv_adapt", "0.25"], "not available in this build yet")
    refuses("a reliability request", ["-reliability", "within"], "not available in this build yet")
    refuses("several -nfeatures values without -cv", ["-nfeatures", "3,4"], "need -cv")
    refuses("-save_r2 without -cv", ["-save_r2"], "need -cv")
    bad = d / "bad_nomask.txt"
    bad.write_text("\n".join(["Subj Run InputFile"] + [f"{s} {r} {f}" for s, g, r, m, f in rows[:6]]) + "\n")
    refuses("datasets with no mask anywhere", [], "need a mask", table=bad)
    short = d / "bad_run.txt"
    short.write_text("\n".join(lines[:-1]) + "\n")                 # drops one run of the last 1D subject
    refuses("a subject missing a run", [], "", table=short)
    before = (d / "out.S.1D").stat().st_mtime_ns
    again = run([prog, "-prefix", out, "-nfeatures", K, "-niter", 5, "-dataTableFile", d / "table.txt", "-progress", "off"],
                expect_ok=False)
    check("refuses: existing outputs without -overwrite", float(not (again.returncode != 0 and "already exists" in again.stdout + again.stderr)), 0.0, "{:.0f}")
    # the refusal comes before the work: it must not have rewritten anything
    check("... and left the earlier output untouched", float((d / "out.S.1D").stat().st_mtime_ns != before), 0.0, "{:.0f}")


# --------------------------------------------------------------------- cross-validation

def np_polar(a):
    u, _, vt = np.linalg.svd(a, full_matrices=False)
    return u @ vt


def np_fisher(r):
    return float(np.tanh(np.arctanh(np.clip(r, -0.9999999, 0.9999999)).mean()))


def np_isc(p, ref):
    return np_fisher(np.array([np.corrcoef(p[c], ref[c])[0, 1] for c in range(p.shape[0])]))


def np_segment_match(a, b, win):
    """BrainIAK-example recipe for one subject against one reference."""
    k, T = a.shape
    nseg = T - win
    za = np.array([a[:, i:i + win].ravel() for i in range(nseg)])
    zb = np.array([b[:, i:i + win].ravel() for i in range(nseg)])
    za = (za - za.mean(1, keepdims=True)); za /= np.maximum(np.linalg.norm(za, axis=1, keepdims=True), 1e-300)
    zb = (zb - zb.mean(1, keepdims=True)); zb /= np.maximum(np.linalg.norm(zb, axis=1, keepdims=True), 1e-300)
    cm = za @ zb.T
    for i in range(nseg):
        for j in range(nseg):
            if abs(i - j) < win and i != j:
                cm[i, j] = -np.inf
    ncand = np.array([1 + sum(1 for j in range(nseg) if abs(i - j) >= win) for i in range(nseg)])
    return float((np.argmax(cm, axis=1) == np.arange(nseg)).mean()), float(np.mean(1.0 / ncand))


def read_scores(path):
    rows = [ln.split("\t") for ln in Path(path).read_text().splitlines() if ln and not ln.startswith("#")]
    head, body = rows[0], rows[1:]
    return head, body


def num(x):
    return float("nan") if x == "NA" else float(x)


def make_cv_world(rng, runs, nsub=8):
    """A world on a common mask with the given run structure."""
    global RUNS
    old = RUNS
    RUNS = runs
    try:
        s_true, subs = make_world(rng, nsub)
        common = np.logical_and.reduce([s["mask"] for s in subs])
        for s in subs:
            sel = common.ravel(order="F")[s["idx"]]
            s["x"], s["idx"] = s["x"][sel], s["idx"][sel]
        return s_true, subs, common
    finally:
        RUNS = old


def case_cv(prog, cli, work, rng, mode, algo, runs, win=5, K_=3, niter=50, label=""):
    global RUNS
    name = f"cv_{mode}_{algo}_{'x'.join(map(str, runs))}{label}"
    print(f"\n[{name}]", flush=True)
    d = work / name
    old = RUNS
    RUNS = runs
    try:
        s_true, subs, common = make_cv_world(rng, runs)
        rows = write_world(d, subs)
        write_table(d / "table.txt", rows, False)
        write_nifti(d / "mask.nii.gz", common.astype(np.float32)[..., None])
        cmd = [prog, "-prefix", d / "out", "-nfeatures", K_, "-algorithm", algo, "-niter", niter, "-tol", 0,
               "-dataTableFile", d / "table.txt", "-mask", d / "mask.nii.gz", "-cv", mode, "-match_window", win,
               "-save_r2", "-progress", "off", "-overwrite"]
        p = run(cmd)
        log = p.stdout + p.stderr
        named = [ln for ln in log.splitlines() if "3dSRM" in ln]
        check("the program name appears on the banner line only", float(len(named) - 1), 0.0, "{:.0f}")
        pre = [preprocess(s)[0] for s in subs]
        n, T = len(pre), pre[0].shape[1]
        folds = []
        c0 = 0
        for L in runs:
            folds.append((c0, L))
            c0 += L
        if len(runs) == 1:                                   # one run: two contiguous blocks
            folds = [(0, T // 2), (T // 2, T - T // 2)]
        head, body = read_scores(d / "out.scores_folds.txt")
        got = {(r[0], r[2]): [num(r[4]), num(r[5]), num(r[6]), num(r[7])] for r in body}
        names = [s["label"] for s in subs]
        fold_names = sorted({r[2] for r in body}, key=lambda x: [r[2] for r in body].index(x))
        want = {}
        vox_ref = [np.zeros(z.shape[0]) for z in pre]
        if mode == "run":
            full = reference_fit(cli, work, name + "_full", pre, K_, niter, algo)
            for fi_, (c0, L) in enumerate(folds):
                tr = [np.delete(z, np.arange(c0, c0 + L), axis=1) for z in pre]
                te = [z[:, c0:c0 + L] for z in pre]
                r = reference_fit(cli, work, f"{name}_f{fi_}", tr, K_, niter, algo)
                P = [r["w"][i].T @ te[i] for i in range(n)]
                for i in range(n):
                    oth = [j for j in range(n) if j != i]
                    ref_u = sum(P[j] for j in oth) / len(oth)
                    if algo == "prob":
                        prec = np.linalg.inv(r["sigma"]) + (1.0 / r["rho2"][oth]).sum() * np.eye(K_)
                        ref_r = np.linalg.solve(prec, sum(P[j] / r["rho2"][j] for j in oth))
                    else:
                        ref_r = ref_u
                    acc, ch = np_segment_match(P[i], ref_u, win)
                    xh = r["w"][i] @ ref_r
                    sse = ((te[i] - xh) ** 2).sum(1)
                    sst = ((te[i] - te[i].mean(1, keepdims=True)) ** 2).sum(1)
                    want[(names[i], fold_names[fi_])] = [np_isc(P[i], ref_u), acc, ch, 1.0 - sse.sum() / sst.sum()]
                    vox_ref[i] += np.where(sst > 0, 1.0 - sse / np.where(sst > 0, sst, 1), 0.0) / len(folds)
                # held-out projections, compared up to a rotation of the gauge, in this fold's columns
                if fi_ == 0:
                    pass
            worst_proj = 0.0
            for i in range(n):
                pj = load_afni(d / f"out.proj.{names[i]}").reshape(K_, T)
                for fi_, (c0, L) in enumerate(folds):
                    tr = [np.delete(z, np.arange(c0, c0 + L), axis=1) for z in pre]
                    r = run_ours(cli, work, f"{name}_f{fi_}", algo, [x.astype(np.float32) for x in tr], K_, niter, False,
                                   extra=["--canon"])
                    Pref = r["w"][i].T @ pre[i][:, c0:c0 + L]
                    worst_proj = max(worst_proj, procrustes_residual(pj[:, c0:c0 + L], Pref))
            check("held-out projections equal the fold models' up to a rotation of the gauge", worst_proj, 2e-3)
        else:
            # leave one subject out: warm start from the full fit without the subject, loso_iter = 5 steps
            full = reference_fit(cli, work, name + "_full", pre, K_, niter, algo)
            for i in range(n):
                oth = [j for j in range(n) if j != i]
                dd = work / f"{name}_drop{i}"
                dd.mkdir(parents=True, exist_ok=True)
                extra = []
                if algo == "prob":
                    full["rho2"][oth].astype(np.float64).tofile(dd / "rho2.f64")
                    full["sigma"].astype(np.float32).tofile(dd / "sigma.f32")
                    extra = ["--init-rho2", str(dd / "rho2.f64"), "--init-sigma", str(dd / "sigma.f32")]
                r = run_ours(cli, work, f"{name}_drop{i}", algo, [pre[j].astype(np.float32) for j in oth], K_, 5, False,
                               init_w=[full["w"][j] for j in oth], extra=extra)
                for fi_, (c0, L) in enumerate(folds):
                    tr = np.delete(pre[i], np.arange(c0, c0 + L), axis=1)
                    te = pre[i][:, c0:c0 + L]
                    s_tr = np.delete(r["s"], np.arange(c0, c0 + L), axis=1)
                    s_te = r["s"][:, c0:c0 + L]
                    W = np_polar(tr @ s_tr.T)
                    Pm = W.T @ te
                    acc, ch = np_segment_match(Pm, s_te, win)
                    sse = ((te - W @ s_te) ** 2).sum(1)
                    sst = ((te - te.mean(1, keepdims=True)) ** 2).sum(1)
                    want[(names[i], fold_names[fi_])] = [np_isc(Pm, s_te), acc, ch, 1.0 - sse.sum() / sst.sum()]
                    vox_ref[i] += np.where(sst > 0, 1.0 - sse / np.where(sst > 0, sst, 1), 0.0) / len(folds)
        errs = np.array([[abs(got[k_][q] - want[k_][q]) for q in range(4)] for k_ in want])
        check("isc, per subject and fold, vs the reference", float(errs[:, 0].max()), 5e-4)
        # the tables carry 7 significant digits, so a fraction like 0.78431372 is rounded in the file
        check("segment-matching accuracy vs the reference", float(errs[:, 1].max()), 1e-6)
        check("... and its chance level", float(errs[:, 2].max()), 1e-6)
        check("held-out R^2 vs the reference", float(errs[:, 3].max()), 5e-4)
        worst_r2map = 0.0
        for i in range(n):
            z, ok, _ = (pre[i], None, None)
            vol = load_afni(d / f"out.r2.{names[i]}").reshape(-1, order="F")
            # the voxels used by this subject, in dataset order
            _, okmask, _ = preprocess(subs[i])
            worst_r2map = max(worst_r2map, float(np.abs(vol[subs[i]["idx"][okmask]] - vox_ref[i]).max()))
        check("voxelwise held-out R^2 maps (-save_r2) vs the reference", worst_r2map, 2e-3)
        # per-subject file: the mean over folds, and the carried covariates, and the factor columns
        head_s, body_s = read_scores(d / "out.scores.txt")
        check("scores.txt: one row per subject with Factor/Level columns", float(not (len(body_s) == len(subs) and head_s[2:4] == ["Factor", "Level"])), 0.0, "{:.0f}")
        mean_isc = np.array([np.mean([want[(names[i], f)][0] for f in fold_names]) for i in range(len(subs))])
        check("scores.txt: the subject mean over folds", float(np.abs(np.array([num(r[5]) for r in body_s]) - mean_isc).max()), 5e-4)
    finally:
        RUNS = old


def case_nfeat_scan(prog, cli, work, rng):
    global RUNS
    print("\n[-nfeatures scan]", flush=True)
    d = work / "scan"
    s_true, subs, common = make_cv_world(rng, RUNS)
    rows = write_world(d, subs)
    write_table(d / "table.txt", rows, False)
    write_nifti(d / "mask.nii.gz", common.astype(np.float32)[..., None])
    base = [prog, "-nfeatures", "2,3,4", "-niter", 40, "-tol", 0, "-dataTableFile", d / "table.txt", "-mask", d / "mask.nii.gz",
            "-cv", "run", "-match_window", 5, "-progress", "off", "-overwrite"]
    p = run(base + ["-prefix", d / "scan"])
    log = p.stdout + p.stderr
    check("the scan reports the best -nfeatures", float("best -nfeatures by held-out isc" not in log), 0.0, "{:.0f}")
    check("the scan writes no model", float((d / "scan.srm").exists() or (d / "scan.S.1D").exists()), 0.0, "{:.0f}")
    head, body = read_scores(d / "scan.nfeat.txt")
    check("nfeat.txt has a row per value", float(len(body) != 3), 0.0, "{:.0f}")
    run([prog, "-prefix", d / "k3", "-nfeatures", 3, "-niter", 40, "-tol", 0, "-dataTableFile", d / "table.txt", "-mask", d / "mask.nii.gz",
         "-cv", "run", "-match_window", 5, "-progress", "off", "-overwrite"])
    _, single = read_scores(d / "k3.scores.txt")
    row3 = [r for r in body if r[0] == "3"][0]
    check("the K=3 row equals the mean over subjects of a single-K run (isc)",
          abs(num(row3[1]) - float(np.mean([num(r[5]) for r in single]))), 1e-6)
    check("... and match", abs(num(row3[3]) - float(np.mean([num(r[6]) for r in single]))), 1e-6)



# --------------------------------------------------------------------- model_in, update_S, apply_W

def table_for(path, rows, subjects, with_mask=False):
    keep = [r for r in rows if r[0] in subjects]
    write_table(path, keep, with_mask)


def case_model_in(prog, cli, work, rng, algo):
    name = f"model_in_{algo}"
    print(f"\n[{name}]", flush=True)
    d = work / name
    runs = [60, 50, 70]
    s_true, subs, common = make_cv_world(rng, runs, nsub=10)
    subs[7]["x"][3, :] = 5.0                                     # a feature that never varies, in a new subject
    old = globals()["RUNS"]
    globals()["RUNS"] = runs
    try:
        rows = write_world(d, subs)
        pre = [preprocess(s)[0] for s in subs]
    finally:
        globals()["RUNS"] = old
    names = [s["label"] for s in subs]
    write_nifti(d / "mask.nii.gz", common.astype(np.float32)[..., None])
    base_names, up1_names, up2_names = names[:6], names[6:8], names[8:]
    table_for(d / "base.txt", rows, base_names)
    table_for(d / "new1.txt", rows, up1_names)
    table_for(d / "new2.txt", rows, up2_names)
    common_args = ["-mask", d / "mask.nii.gz", "-progress", "off", "-overwrite"]
    K_, T = 3, sum(runs)

    run([prog, "-prefix", d / "base", "-nfeatures", K_, "-algorithm", algo, "-niter", 60, "-tol", 0, "-save_W",
         "-dataTableFile", d / "base.txt"] + common_args)
    hdr, subj, _, blocks = parse_srm_file(d / "base.srm")
    S = blocks["S"]
    base_bytes = (d / "base.srm").read_bytes()

    # -- A: a W for each new subject against the fixed S
    p = run([prog, "-prefix", d / "mi", "-model_in", d / "base.srm", "-dataTableFile", d / "new1.txt"] + common_args)
    log = p.stdout + p.stderr
    check("A: the program name appears on the banner line only", float(len([l for l in log.splitlines() if "3dSRM" in l]) - 1), 0.0, "{:.0f}")
    check("A: no new model and no S are written", float((d / "mi.srm").exists() or (d / "mi.S.1D").exists()), 0.0, "{:.0f}")
    check("A: the model file is untouched", float((d / "base.srm").read_bytes() != base_bytes), 0.0, "{:.0f}")
    worst, worst_r2 = 0.0, 0.0
    rho2 = None
    if algo == "prob":
        rho2 = np.array([float(x) for x in Path(d / "mi.rho2.1D").read_text().splitlines()[1:]])
    for j, nm in enumerate(up1_names):
        z = pre[names.index(nm)]
        W = np_polar(z @ S.T)
        P = load_afni(d / f"mi.proj.{nm}").reshape(K_, T)
        worst = max(worst, float(np.abs(P - W.T @ z).max()))
        if rho2 is not None:
            r2 = ((z - W @ S) ** 2).mean()
            worst_r2 = max(worst_r2, abs(rho2[j] - r2) / r2)
    check("A: new subjects' data in the shared space vs numpy (W = polar(X S'))", worst, 2e-4)
    if rho2 is not None:
        check("A: their noise variances vs numpy", worst_r2, 1e-5)

    # -- B: -cv run scores a new subject on runs neither its W nor S saw
    run([prog, "-prefix", d / "cvmi", "-model_in", d / "base.srm", "-cv", "run", "-match_window", 5, "-save_r2",
         "-dataTableFile", d / "new1.txt"] + common_args)
    head, body = read_scores(d / "cvmi.scores_folds.txt")
    got = {(r[0], r[2]): [num(r[4]), num(r[5]), num(r[6]), num(r[7])] for r in body}
    fold_names = sorted({r[2] for r in body})
    folds, c0 = [], 0
    for L in runs:
        folds.append((c0, L))
        c0 += L
    errs = []
    for nm in up1_names:
        z = pre[names.index(nm)]
        for fi_, (c0, L) in enumerate(folds):
            tr = np.delete(z, np.arange(c0, c0 + L), axis=1)
            te = z[:, c0:c0 + L]
            s_tr = np.delete(S, np.arange(c0, c0 + L), axis=1)
            s_te = S[:, c0:c0 + L]
            W = np_polar(tr @ s_tr.T)
            Pm = W.T @ te
            acc, ch = np_segment_match(Pm, s_te, 5)
            sse = ((te - W @ s_te) ** 2).sum(1)
            sst = ((te - te.mean(1, keepdims=True)) ** 2).sum(1)
            wnt = [np_isc(Pm, s_te), acc, ch, 1.0 - sse.sum() / sst.sum()]
            errs.append([abs(got[(nm, fold_names[fi_])][q] - wnt[q]) for q in range(4)])
    errs = np.array(errs)
    check("B: held-out isc vs numpy", float(errs[:, 0].max()), 5e-4)
    check("B: segment matching vs numpy", float(errs[:, 1].max()), 1e-6)
    check("B: ... and its chance level", float(errs[:, 2].max()), 1e-6)
    check("B: held-out R^2 vs numpy", float(errs[:, 3].max()), 5e-4)

    # -- C: -update_S against the library driver, and the chain of two updates
    def ref_update(model_file, tag, upd_names, niter):
        rows_ = [pre[names.index(n)] for n in upd_names]
        files = []
        for n, z in zip(upd_names, rows_):
            f = d / f"{tag}_{n}.f32"
            np.ascontiguousarray(z, dtype=np.float32).tofile(f)
            files.append(f)
        out = d / f"{tag}_ref"
        run_load_update(cli, model_file, out, files, [z.shape[0] for z in rows_], niter)
        S_ref = np.fromfile(f"{out}.u_S.f32", dtype=np.float32).reshape(K_, T).astype(np.float64)
        a_ref = np.fromfile(f"{out}.u_asum.f64").reshape(K_, T)
        info = open(f"{out}.u_info.txt").read().split()
        return S_ref, a_ref, float(info[2]), int(info[3]), int(info[4])

    p = run([prog, "-prefix", d / "up1", "-model_in", d / "base.srm", "-update_S", "-niter", 80, "-tol", 0,
             "-cv", "run", "-match_window", 5, "-dataTableFile", d / "new1.txt"] + common_args)
    log = p.stdout + p.stderr
    S_ref, a_ref, drift_ref, nupd_ref, neff_ref = ref_update(d / "base.srm", "u1", up1_names, 80)
    S_got = np.loadtxt(d / "up1.S.1D", comments="#").T
    check("C: updated S vs the library driver (relative)", float(np.abs(S_got - S_ref).max() / np.abs(S_ref).max()), 2e-4)
    h1, s1, _, b1 = parse_srm_file(d / "up1.srm")
    check("C: updated model: S in the file vs the library driver (relative)",
          float(np.abs(b1["S"] - S_ref).max() / np.abs(S_ref).max()), 2e-4)
    check("C: updated model: a_sum vs the library driver (relative)",
          float(np.abs(b1["a_sum"] - a_ref).max() / np.abs(a_ref).max()), 2e-4)
    check("C: updated model: subjects and update count in the header",
          float(not (len(s1) == 8 and int(h1["nsubj"][0]) == 8 and int(h1["nupdate"][0]) == nupd_ref
                     and int(h1["nsubj_eff"][0]) == neff_ref and [x[0] for x in s1] == names[:8])), 0.0, "{:.0f}")
    check("C: updated model: drift vs the library driver", abs(float(h1["drift"][0]) - drift_ref) / max(drift_ref, 1e-12), 2e-2)
    # scores are computed from the model as it was, before the update
    cvmi = (d / "cvmi.scores.txt").read_text()
    check("C: -cv scores with -update_S equal those without (computed before the update)",
          float((d / "up1.scores.txt").read_text() != cvmi), 0.0, "{:.0f}")

    run([prog, "-prefix", d / "up2", "-model_in", d / "up1.srm", "-update_S", "-niter", 80, "-tol", 0,
         "-dataTableFile", d / "new2.txt"] + common_args)
    S_ref2, a_ref2, drift_ref2, nupd_ref2, neff_ref2 = ref_update(d / "up1.srm", "u2", up2_names, 80)
    h2, s2, _, b2 = parse_srm_file(d / "up2.srm")
    check("C: second update: S vs the library driver (relative)", float(np.abs(b2["S"] - S_ref2).max() / np.abs(S_ref2).max()), 2e-4)
    check("C: second update: 10 subjects, 2 updates",
          float(not (len(s2) == 10 and int(h2["nupdate"][0]) == 2 and nupd_ref2 == 2)), 0.0, "{:.0f}")

    # -- errors and protection
    def fails(args, needle, label):
        pr = run([prog, "-progress", "off"] + args, expect_ok=False)
        out = pr.stdout + pr.stderr
        check(label, float(pr.returncode == 0 or needle not in out), 0.0, "{:.0f}")
        if pr.returncode == 0 or needle not in out:
            print(out)
    fails(["-prefix", d / "e1", "-model_in", d / "base.srm", "-update_S", "-dataTableFile", d / "base.txt", "-mask", d / "mask.nii.gz"],
          "already one of the", "E: -update_S refuses a subject that is already in the model")
    fails(["-prefix", d / "e2", "-model_in", d / "base.srm", "-nfeatures", 3, "-dataTableFile", d / "new1.txt", "-mask", d / "mask.nii.gz"],
          "-nfeatures is defined by the model", "E: -nfeatures with -model_in is refused")
    fails(["-prefix", d / "e3", "-model_in", d / "base.srm", "-cv", "loso", "-dataTableFile", d / "new1.txt", "-mask", d / "mask.nii.gz"],
          "-cv loso", "E: -cv loso with -model_in is refused")
    fails(["-prefix", d / "e4", "-update_S", "-nfeatures", 3, "-dataTableFile", d / "new1.txt", "-mask", d / "mask.nii.gz"],
          "-update_S needs -model_in", "E: -update_S alone is refused")
    fails(["-prefix", d / "e5", "-model_in", d / "base.srm", "-niter", 5, "-dataTableFile", d / "new1.txt", "-mask", d / "mask.nii.gz"],
          "apply only to -update_S", "E: -niter without -update_S is refused")
    short = [r for r in rows if r[0] in up1_names and r[2] != "r3"]
    write_table(d / "short.txt", short, False)
    fails(["-prefix", d / "e6", "-model_in", d / "base.srm", "-dataTableFile", d / "short.txt", "-mask", d / "mask.nii.gz"],
          "do not fit the model", "E: new data with a different run structure are refused with one clear message")
    no_ow = ["-prefix", d / "base", "-model_in", d / "base.srm", "-update_S", "-niter", 20,
             "-dataTableFile", d / "new1.txt", "-mask", d / "mask.nii.gz"]
    fails(no_ow, "already exists", "E: -update_S does not overwrite the model it started from without -overwrite")
    check("E: ... and the model file is still intact", float((d / "base.srm").read_bytes() != base_bytes), 0.0, "{:.0f}")
    run([prog, "-progress", "off"] + no_ow + ["-overwrite"])
    h3, s3, _, _ = parse_srm_file(d / "base.srm")
    check("E: with -overwrite the same name is updated in place", float(len(s3) != 8), 0.0, "{:.0f}")


def case_apply_W(prog, cli, work, rng, zs):
    name = f"apply_W_{zs}"
    print(f"\n[{name}]", flush=True)
    d = work / name
    runs = [60, 50, 70]
    s_true, subs, common = make_cv_world(rng, runs, nsub=5)
    old = globals()["RUNS"]
    globals()["RUNS"] = runs
    try:
        rows = write_world(d, subs)
    finally:
        globals()["RUNS"] = old
    names = [s["label"] for s in subs]
    write_table(d / "base.txt", rows, False)
    write_nifti(d / "mask.nii.gz", common.astype(np.float32)[..., None])
    K_ = 3
    run([prog, "-prefix", d / "base", "-nfeatures", K_, "-niter", 40, "-tol", 0, "-save_W", "-zscore", zs, "-mask", d / "mask.nii.gz",
         "-dataTableFile", d / "base.txt", "-progress", "off", "-overwrite"])
    # other data for the same people: different lengths, even between subjects
    task_runs = {nm: ([40, 55] if i % 2 == 0 else [45, 50]) for i, nm in enumerate(names)}
    task, trows = {}, []
    for i, (s, nm) in enumerate(zip(subs, names)):
        v = len(s["idx"])
        x = rng.standard_normal((v, sum(task_runs[nm]))) * 3.0 + 100.0
        task[nm] = x
        c0 = 0
        mask_file = d / "mask.nii.gz"
        for r, L in enumerate(task_runs[nm]):
            vol = np.zeros((NX, NY, NZ, L), dtype=np.float32)
            flat = vol.reshape(-1, L, order="F")
            flat[s["idx"], :] = x[:, c0:c0 + L]
            write_nifti(d / f"{nm}_task{r + 1}.nii.gz", flat.reshape((NX, NY, NZ, L), order="F"))
            trows.append((nm, s["group"], f"r{r + 1}", mask_file, d / f"{nm}_task{r + 1}.nii.gz"))
            c0 += L
    write_table(d / "task.txt", trows, False)
    p = run([prog, "-prefix", d / "task", "-apply_W", d / "base", "-dataTableFile", d / "task.txt", "-progress", "off", "-overwrite"])
    log = p.stdout + p.stderr
    check("the program name appears on the banner line only", float(len([l for l in log.splitlines() if "3dSRM" in l]) - 1), 0.0, "{:.0f}")
    worst = 0.0
    for s, nm in zip(subs, names):
        Wvol = load_afni(d / f"base.W.{nm}")                      # NX x NY x NZ x K
        Wm = Wvol.reshape(-1, K_, order="F")[s["idx"]] if Wvol.ndim == 4 else None
        Wm = np.stack([Wvol[..., c].reshape(-1, order="F")[s["idx"]] for c in range(K_)], axis=1)
        x = task[nm].astype(np.float32).astype(np.float64)
        if zs == "run":
            z, c0 = np.zeros_like(x), 0
            for L in task_runs[nm]:
                seg = x[:, c0:c0 + L]
                z[:, c0:c0 + L] = (seg - seg.mean(1, keepdims=True)) / seg.std(1, keepdims=True)
                c0 += L
        else:
            z = x - x.mean(1, keepdims=True)
        Tn = x.shape[1]
        P = load_afni(d / f"task.proj.{nm}").reshape(K_, Tn)
        worst = max(worst, float(np.abs(P - Wm.T @ z).max()))
    check("projected new data vs numpy (own W, z-scored as the model's were)", worst, 2e-4)
    pr = run([prog, "-prefix", d / "bad", "-apply_W", d / "base", "-dataTableFile", d / "task.txt", "-mask", d / "mask.nii.gz",
              "-progress", "off", "-overwrite"], expect_ok=False)
    check("-mask with -apply_W is refused", float(pr.returncode == 0), 0.0, "{:.0f}")
    gone = [r for r in trows if r[0] != names[0]]
    gone = [("stranger" if r[0] == names[1] else r[0],) + tuple(r[1:]) for r in gone]
    write_table(d / "stranger.txt", gone, False)
    pr = run([prog, "-prefix", d / "bad2", "-apply_W", d / "base", "-dataTableFile", d / "stranger.txt", "-progress", "off", "-overwrite"], expect_ok=False)
    check("a subject missing from the model is refused", float(pr.returncode == 0 or "not in the model" not in pr.stdout + pr.stderr), 0.0, "{:.0f}")
    pr = run([prog, "-prefix", d / "task", "-apply_W", d / "base", "-dataTableFile", d / "task.txt", "-progress", "off"], expect_ok=False)
    check("existing outputs are not overwritten", float(pr.returncode == 0), 0.0, "{:.0f}")


def case_1D_apply_W(prog, cli, work, rng):
    print("\n[1D subjects: dropped features and -apply_W]", flush=True)
    d = work / "apply_W_1D"
    d.mkdir(parents=True, exist_ok=True)
    T1, K_ = 90, 3
    s = rng.standard_normal((K_, T1))
    rows = []
    cols = {}
    for i in range(4):
        v = 12 + 3 * i
        x = (np.linalg.qr(rng.standard_normal((v, K_)))[0] @ s + 0.5 * rng.standard_normal((v, T1))) * 2.0 + 10.0
        if i == 1:
            x[4, :] = 3.0                                        # a column that never varies
        cols[f"s{i + 1}"] = x
        np.savetxt(d / f"s{i + 1}.1D", x.T, fmt="%.9g")
        rows.append(f"s{i + 1} g {d}/s{i + 1}.1D")
    (d / "base.txt").write_text("Subj Group InputFile\n" + "\n".join(rows) + "\n")
    run([prog, "-prefix", d / "base", "-nfeatures", K_, "-niter", 30, "-tol", 0, "-save_W", "-proj_format", "1D",
         "-dataTableFile", d / "base.txt", "-progress", "off", "-overwrite"])
    W2 = np.loadtxt(d / "base.W.s2.1D", comments="#")
    check("W of a 1D subject has a row for every column of the input", float(W2.shape[0] != cols["s2"].shape[0]), 0.0, "{:.0f}")
    check("... and the dropped column's row is zero", float(np.abs(W2[4]).max()), 0.0)
    new = {k_: rng.standard_normal((v.shape[0], 70)) + 5.0 for k_, v in cols.items()}
    rows = []
    for k_, x in new.items():
        np.savetxt(d / f"{k_}_new.1D", x.T, fmt="%.9g")
        rows.append(f"{k_} g {d}/{k_}_new.1D")
    (d / "new.txt").write_text("Subj Group InputFile\n" + "\n".join(rows) + "\n")
    run([prog, "-prefix", d / "new", "-apply_W", d / "base", "-proj_format", "1D", "-dataTableFile", d / "new.txt",
         "-progress", "off", "-overwrite"])
    worst = 0.0
    for k_, x in new.items():
        W = np.loadtxt(d / f"base.W.{k_}.1D", comments="#")
        xs = x.astype(np.float32).astype(np.float64)
        z = (xs - xs.mean(1, keepdims=True)) / xs.std(1, keepdims=True)
        P = np.loadtxt(d / f"new.proj.{k_}.1D", comments="#").T
        worst = max(worst, float(np.abs(P - W.T @ z).max()))
    check("1D: projected new data vs numpy", worst, 2e-4)


# --------------------------------------------------------------------- -axis item

ITEM_RUNS = ["r1", "r2", "r3"]
ITEM_CATS = ["face", "house", "tool"]


def afni_labels(prefix):
    """Volume labels of an AFNI +orig dataset, read from the .HEAD text (BRICK_LABS)."""
    head = Path(str(prefix) + "+orig.HEAD").read_text()
    i = head.index("name = BRICK_LABS")
    m = re.search(r"count = (\d+)\s*\n'", head[i:])
    text = head[i + m.end():i + m.end() + int(m.group(1)) - 1]
    return text.split("~")[:-1] if text.endswith("~") else text.split("~")


def expected_keys(table_rows, repeat_col=False):
    """The column keys and their order, from the help text: subjects and runs in order of first appearance;
    each subject's rows in run order, then table order; an item a subject has more than once is keyed
    Item:1..n in that order (Item:R with a Repeat column); the first subject's order is the model's."""
    runs = []
    for r in table_rows:
        if r["run"] is not None and r["run"] not in runs:
            runs.append(r["run"])
    subs = []
    for r in table_rows:
        if r["subj"] not in subs:
            subs.append(r["subj"])
    keyed = {}
    for s in subs:
        mine = [(runs.index(r["run"]) if r["run"] is not None else 0, pos, r) for pos, r in enumerate(table_rows) if r["subj"] == s]
        mine.sort(key=lambda t: (t[0], t[1]))
        total = {}
        for _, _, r in mine:
            total[r["item"]] = total.get(r["item"], 0) + 1
        seen, out = {}, []
        for _, _, r in mine:
            seen[r["item"]] = seen.get(r["item"], 0) + 1
            if total[r["item"]] > 1:
                out.append((f"{r['item']}:{r['rep'] if repeat_col else seen[r['item']]}", r))
            else:
                out.append((r["item"], r))
        keyed[s] = out
    return subs, [k for k, _ in keyed[subs[0]]], keyed


def make_item_world(rng, nitem=24, nrep=6, nsub=6, n1d=2, use_runs=True, k=K, shared_runs=False, same_run=False,
                    run_offsets=True):
    """Items with a latent response each; a repeat shares its item's response. Each subject has its own random
    run assignment (or, with shared_runs, everyone has the same one), and, with run_offsets, a per-voxel offset
    in every run (betas of per-run GLMs; without them, betas of one GLM).
    With same_run, the repeats of it01 fall in the same run. Repeats carry labels 'a', 'b' in reverse run order."""
    s_item = rng.standard_normal((k, nitem)) * np.linspace(1.6, 0.7, k)[:, None]
    names = [f"it{j + 1:02d}" for j in range(nitem)]

    def assignment():
        occ = []
        for j in range(nitem):
            if not use_runs:
                occ.append((j, None, None))
                continue
            rr = rng.permutation(len(ITEM_RUNS))
            if j < nrep:
                r1, r2 = (int(rr[0]), int(rr[0])) if (same_run and j == 0) else (int(rr[0]), int(rr[1]))
                first, second = sorted([r1, r2])
                occ += [(j, first, "b"), (j, second, "a")]
            else:
                occ.append((j, int(rr[0]), None))
        return occ
    common = assignment() if shared_runs else None
    subs = []
    for i in range(nsub):
        occ = list(common) if shared_runs else assignment()
        if i >= nsub - n1d:
            v = 12 + 3 * i
            mask, idx = None, None
        else:
            mask = np.zeros((NX, NY, NZ), dtype=bool)
            if i % 2 == 0:
                mask[1:9, 1:9, 1:5] = True               # 256 voxels
            else:
                mask[2:8, 1:9, 1:5] = True               # 192: masks differ between subjects
            idx = np.flatnonzero(mask.ravel(order="F"))
            v = len(idx)
        q, _ = np.linalg.qr(rng.standard_normal((v, k)))
        scale = rng.uniform(0.5, 3.0, (v, 1))
        offs = {rr: rng.uniform(-40, 40, (v, 1)) for rr in range(len(ITEM_RUNS))}
        cols = []
        for j, rr, _ in occ:
            x = q @ s_item[:, j:j + 1] + 0.5 * rng.standard_normal((v, 1))
            cols.append(x * scale + 200 + (offs[rr] if (rr is not None and run_offsets) else 0.0))
        subs.append({"label": f"s{i + 1:02d}", "group": "ctl" if i % 2 == 0 else "pat", "mask": mask, "idx": idx,
                     "occ": occ, "vals": cols, "is1D": i >= nsub - n1d})
    if use_runs:                                        # flat in one run for one subject: dropped under -zscore run
        s0 = subs[0]
        for c, (j, rr, _) in enumerate(s0["occ"]):
            if rr == 1:
                s0["vals"][c][4] = 9.0
    return names, subs


def write_item_world(d, names, subs, rng, layout_cycle=("per_run", "one_file", "per_item")):
    """Write the inputs (datasets read through [n] selectors in three layouts, or one-row .1D files) and
    return the table rows, shuffled."""
    d.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, s in enumerate(subs):
        lab = s["label"]
        occ_rows = []
        if s["is1D"]:
            for c in range(len(s["occ"])):
                f = d / f"{lab}_c{c:03d}.1D"
                np.savetxt(f, s["vals"][c].T, fmt="%.8g")                 # one row: the item's features
                occ_rows.append((c, str(f)))
        else:
            mask_file = d / f"{lab}_mask.nii.gz"
            write_nifti(mask_file, s["mask"].astype(np.float32)[..., None])
            layout = layout_cycle[i % len(layout_cycle)]
            groups = {}
            for c, (j, rr, _) in enumerate(s["occ"]):
                key = "all" if layout == "one_file" else (f"c{c:03d}" if layout == "per_item" else f"run{rr}")
                groups.setdefault(key, []).append(c)
            for key, cs in groups.items():
                order = list(rng.permutation(cs))                         # volumes in a scrambled order
                vol = np.zeros((NX, NY, NZ, len(order)), dtype=np.float32)
                flat = vol.reshape(-1, len(order), order="F")
                for b, c in enumerate(order):
                    flat[s["idx"], b] = s["vals"][c][:, 0]
                vol = flat.reshape((NX, NY, NZ, len(order)), order="F")
                f = d / f"{lab}_{key}.nii.gz"
                write_nifti(f, vol)
                for b, c in enumerate(order):
                    occ_rows.append((c, f"{f}[{b}]" if len(order) > 1 else str(f)))
        for c, path in occ_rows:
            j, rr, rep = s["occ"][c]
            rows.append({"subj": lab, "group": s["group"], "item": names[j], "c": c, "rep": rep or "-",
                         "run": ITEM_RUNS[rr] if rr is not None else None, "cat": ITEM_CATS[j % len(ITEM_CATS)],
                         "mask": "-" if s["is1D"] else str(d / f"{lab}_mask.nii.gz"), "file": path})
    order = rng.permutation(len(rows))
    return [rows[p] for p in order]


def write_item_table(path, rows, with_run=True, with_repeat=False, with_cat=False):
    head = "Subj Group Item " + ("Run " if with_run else "") + ("Repeat " if with_repeat else "") + \
           ("Category " if with_cat else "") + "Mask InputFile"
    lines = [head]
    for r in rows:
        lines.append(f"{r['subj']} {r['group']} {r['item']} " + (f"{r['run']} " if with_run else "") +
                     (f"{r['rep']} " if with_repeat else "") + (f"{r['cat']} " if with_cat else "") + f"{r['mask']} {r['file']}")
    Path(path).write_text("\n".join(lines) + "\n")


def item_preprocess(sub, keyed_rows, keys, zmode="all"):
    """The subject's features x keys matrix: z-scored over all items ('all'), over the items of each of the
    subject's own runs ('run'), or as given ('none'); features that never vary (in a run, for 'run') dropped."""
    x = np.zeros((sub["vals"][0].shape[0], len(keys)))
    run_of = np.zeros(len(keys), dtype=int)
    for key, r in keyed_rows:
        j = keys.index(key)
        x[:, j] = sub["vals"][r["c"]][:, 0].astype(np.float32)
        run_of[j] = ITEM_RUNS.index(r["run"]) if (zmode == "run" and r["run"] is not None) else 0
    z = np.zeros_like(x)
    bad = np.zeros(x.shape[0], dtype=bool)
    for rr in np.unique(run_of):
        cols = run_of == rr
        seg = x[:, cols]
        sd = seg.std(axis=1)
        bad |= ~(sd > 1e-7 * np.abs(seg.mean(axis=1)))
        z[:, cols] = seg if zmode == "none" else (seg - seg.mean(axis=1, keepdims=True)) / np.where(sd > 0, sd, 1.0)[:, None]
    return z[~bad], ~bad


def s1d_keys(path):
    for ln in Path(path).read_text().splitlines():
        if ln.startswith("# items:"):
            return ln.split()[2:]
    return None


def item_fails(prog, d, args, needle, label, table):
    pr = run([prog, "-progress", "off", "-prefix", d / "bad", "-overwrite", "-nfeatures", K, "-dataTableFile", table] + args,
             expect_ok=False)
    out = pr.stdout + pr.stderr
    check(label, float(pr.returncode == 0 or needle not in out), 0.0, "{:.0f}")
    if pr.returncode == 0 or needle not in out:
        print(out)


def case_items(prog, cli, work, rng, algo, use_runs=True, zmode="all"):
    name = f"items_{algo}_{'runs' if use_runs else 'noruns'}_z{zmode}"
    print(f"\n[{name}]", flush=True)
    d = work / name
    names, subs = make_item_world(rng, use_runs=use_runs, same_run=use_runs)
    rows = write_item_world(d, names, subs, rng)
    write_item_table(d / "table.txt", rows, with_run=use_runs)
    niter = 60
    zopt = [] if zmode == "all" else ["-zscore", zmode]
    p = run([prog, "-prefix", d / "out", "-axis", "item", "-nfeatures", K, "-algorithm", algo, "-niter", niter, "-tol", 0,
             "-dataTableFile", d / "table.txt", "-progress", "line", "-overwrite"] + zopt)
    log = p.stdout + p.stderr
    check("the program name appears on the banner line only", float(len([l for l in log.splitlines() if "3dSRM" in l]) - 1), 0.0, "{:.0f}")
    check("the reading progress counts items", float("reading:" not in log or "items" not in log), 0.0, "{:.0f}")

    sub_order, keys, keyed = expected_keys(rows)
    by_label = {s["label"]: s for s in subs}
    pre = [item_preprocess(by_label[lab], keyed[lab], keys, zmode) for lab in sub_order]
    ref = run_ours(cli, work, name + "_ref", algo, [z.astype(np.float32) for z, _ in pre], K, niter, zmode == "none",
                   extra=["--canon"])

    S = np.loadtxt(d / "out.S.1D").T
    check("item keys in S.1D, in the order the help describes (repeats in one run by row)",
          float(s1d_keys(d / "out.S.1D") != keys), 0.0, "{:.0f}")
    check(f"S.1D (one row per item, -zscore {zmode}) vs the library driver", relfro(S, ref["s"]), 2e-4)
    rows_fit = [ln.split("\t") for ln in (d / "out.fit.txt").read_text().splitlines() if ln and not ln.startswith(("#", "Subj"))]
    check("subjects in order of first appearance", float([r[0] for r in rows_fit] != sub_order), 0.0, "{:.0f}")
    check("features used per subject (constant features dropped)",
          float(sum(abs(int(r[3]) - ok.sum()) for r, (_, ok) in zip(rows_fit, pre))), 0.0, "{:.0f}")
    worst, labs_ok = 0.0, True
    for i, lab in enumerate(sub_order):
        z, _ = pre[i]
        if by_label[lab]["is1D"]:
            continue
        proj = load_afni(d / f"out.proj.{lab}").reshape(K, len(keys))
        mu = z.mean(1, keepdims=True) if zmode == "none" else 0.0
        worst = max(worst, relfro(proj, ref["w"][i].T @ (z - mu)))
        labs_ok &= afni_labels(d / f"out.proj.{lab}") == keys
    check("projected items vs W'X from the driver, worst subject", worst, 5e-3)
    check("projected datasets are labelled with the item keys", float(not labs_ok), 0.0, "{:.0f}")
    hdr, subj, items, blocks = parse_srm_file(d / "out.srm")
    check("model file: axis item, -zscore, and the keys in column order",
          float(hdr.get("axis", [""])[0] != "item" or hdr["zscore"][0] != zmode or items != keys or int(hdr["nitem"][0]) != len(keys)),
          0.0, "{:.0f}")
    if not use_runs or zmode != "all":
        return

    # the row order of the table matters only through the order in which subjects first appear (which
    # sets each one's random start) and runs first appear (which orders repeats): keep the shortest
    # prefix that fixes both, and permute the rest so that no two rows of the same item and subject
    # trade places (that would renumber repeats in one run)
    seen_s, seen_r, m = set(), set(), 0
    while len(seen_s) < len(sub_order) or len(seen_r) < len(ITEM_RUNS):
        seen_s.add(rows[m]["subj"]); seen_r.add(rows[m]["run"]); m += 1
    rest = [rows[m + p] for p in rng.permutation(len(rows) - m)]
    same_run_pairs = [r for r in rows if r["item"] == "it01"]
    rest = [r for r in rest if r["item"] != "it01"]
    rows2 = rows[:m] + [r for r in same_run_pairs if not any(r is q for q in rows[:m])] + rest
    write_item_table(d / "table2.txt", rows2)
    run([prog, "-prefix", d / "out2", "-axis", "item", "-nfeatures", K, "-algorithm", algo, "-niter", niter, "-tol", 0,
         "-dataTableFile", d / "table2.txt", "-progress", "off", "-overwrite"])
    keys2 = s1d_keys(d / "out2.S.1D")
    S2 = np.loadtxt(d / "out2.S.1D").T
    perm = [keys2.index(k_) for k_ in keys]
    check("the permuted table gives the columns in another order", float(keys2 == keys), 0.0, "{:.0f}")
    check("table row order does not matter (S matched by item key)", relfro(S2[:, perm], S), 2e-5)

    # a Repeat column names repeats: here 'b' is the earlier run's, so the keys follow the labels, not the runs
    write_item_table(d / "table_rep.txt", rows, with_repeat=True)
    run([prog, "-prefix", d / "rep", "-axis", "item", "-nfeatures", K, "-niter", 5, "-dataTableFile", d / "table_rep.txt",
         "-progress", "off", "-overwrite"])
    _, keys_rep, _ = expected_keys(rows, repeat_col=True)
    check("a Repeat column keys repeats as Item:R", float(s1d_keys(d / "rep.S.1D") != keys_rep or "it02:a" not in keys_rep), 0.0, "{:.0f}")

    # repeats without a Run column are allowed, numbered by table order
    write_item_table(d / "table_norun.txt", rows, with_run=False)
    run([prog, "-prefix", d / "norun", "-axis", "item", "-nfeatures", K, "-niter", 5, "-dataTableFile", d / "table_norun.txt",
         "-progress", "off", "-overwrite"])
    _, keys_nr, _ = expected_keys([dict(r, run=None) for r in rows])
    check("repeats without a Run column: numbered by table order", float(s1d_keys(d / "norun.S.1D") != keys_nr), 0.0, "{:.0f}")

    # refusals, each naming what to do
    def table_with(rows_, fname, **kw):
        write_item_table(d / fname, rows_, **kw)
        return d / fname
    T0 = d / "table.txt"
    ax = ["-axis", "item"]
    for opt, needle in ((["-censor", T0], "-censor drops time points"), (["-zcensor"], "-zcensor is for time"),
                        (["-cv_blocks", 3], "items have no order"), (["-match_window", 5], "identified on its own"),
                        (["-time_labels", "x", T0], "table column"), (["-repeated_stimulus"], "recur across runs"),
                        (["-run_coupling", 10], "-axis time only")):
        item_fails(prog, d, ax + opt, needle, f"E: {opt[0]} with -axis item is refused, saying why", T0)
    item_fails(prog, d, ["-item_column", "Item"], "need -axis item", "E: -item_column without -axis item", T0)
    item_fails(prog, d, ["-fold_column", "Run", "-cv", "run"], "-fold_column is for -axis item", "E: -fold_column with -axis time", T0)
    item_fails(prog, d, ax + ["-fold_column", "Run"], "needs -cv", "E: -fold_column without -cv", T0)
    item_fails(prog, d, ax + ["-cv", "run"], "so the runs cannot be folds", "E: randomized runs cannot be -cv folds, saying what to do", T0)
    item_fails(prog, d, ax + ["-zscore", "run"], "needs at least 3", "E: -zscore run with a run of fewer than 3 items",
               table_with([dict(r, run=("r9" if (r["subj"] == sub_order[1] and r["item"] in ("it23", "it24")) else r["run"])) for r in rows], "tiny.txt"))
    item_fails(prog, d, ax + ["-zscore", "run"], "no 'Run' column", "E: -zscore run without a Run column", table_with(rows, "norun2.txt", with_run=False))
    s2 = sub_order[1]
    miss = [r for r in rows if not (r["subj"] == s2 and r["item"] == "it20")]
    item_fails(prog, d, ax, f"Subj {s2} is missing Item it20", "E: a missing item names the Subj and Item", table_with(miss, "miss.txt"))
    extra = rows + [dict(next(r for r in rows if r["subj"] == s2 and r["item"] == "it20"), item="it99")]
    item_fails(prog, d, ax, f"Subj {s2} has Item it99, which Subj", "E: an extra item names the Subj and Item", table_with(extra, "extra.txt"))
    one = next(r for r in rows if r["subj"] == s2 and r["item"] == "it02")
    fewer = [r for r in rows if r is not one]
    item_fails(prog, d, ax, f"Subj {s2} has Item it02 1 time", "E: a repeat missing for one subject gives both counts", table_with(fewer, "fewer.txt"))
    duprep = [dict(r, rep="a") if (r["subj"] == s2 and r["item"] == "it02") else r for r in rows]
    item_fails(prog, d, ax, "with Repeat a twice", "E: a Repeat label used twice for one item", table_with(duprep, "duprep.txt", with_repeat=True))
    colon = [dict(r, item="a:b") if r is rows[0] else r for r in rows]
    item_fails(prog, d, ax, "may not contain ':'", "E: an item label with ':' is refused", table_with(colon, "colon.txt"))
    ds_row = next(r for r in rows if r["file"].endswith("]"))
    multi = [dict(r, file=r["file"].rsplit("[", 1)[0] + "[0,1]") if r is ds_row else r for r in rows]
    item_fails(prog, d, ax, "selects 2 volumes", "E: a selector picking two volumes is refused", table_with(multi, "multi.txt"))
    zero_sub = next(s for s in subs if not s["is1D"])
    zf = d / "zero_item.nii.gz"
    write_nifti(zf, np.zeros((NX, NY, NZ, 1), dtype=np.float32))
    zrow = next(r for r in rows if r["subj"] == zero_sub["label"])
    zero = [dict(r, file=str(zf)) if r is zrow else r for r in rows]
    item_fails(prog, d, ax, "is all zero in the mask", "E: an all-zero item names the Subj and item", table_with(zero, "zero.txt"))
    one_d = next(s for s in subs if s["is1D"])
    two = d / "two_rows.1D"
    np.savetxt(two, np.vstack([one_d["vals"][0].T, one_d["vals"][1].T]), fmt="%.8g")
    r1d = next(r for r in rows if r["subj"] == one_d["label"])
    tworow = [dict(r, file=str(two)) if r is r1d else r for r in rows]
    item_fails(prog, d, ax, "a .1D input is one row", "E: a .1D item with two rows is refused", table_with(tworow, "tworow.txt"))


def np_item_identify(a, b):
    """SRM_item_identify written from its contract: item j of a goes to the item of b it correlates with best
    across components (ties to the lower index); chance 1/n."""
    za = a - a.mean(0, keepdims=True); za /= np.where(np.linalg.norm(za, axis=0) > 0, np.linalg.norm(za, axis=0), 1)
    zb = b - b.mean(0, keepdims=True); zb /= np.where(np.linalg.norm(zb, axis=0) > 0, np.linalg.norm(zb, axis=0), 1)
    cm = za.T @ zb
    return float((np.argmax(cm, axis=1) == np.arange(a.shape[1])).mean()), 1.0 / a.shape[1]


def case_items_cv(prog, cli, work, rng, mode, algo, folds_by, k=5, niter=50):
    """-cv on items: folds are the runs (shared by everyone here) or the levels of a category column."""
    name = f"items_cv_{mode}_{algo}_{folds_by}"
    print(f"\n[{name}]", flush=True)
    d = work / name
    names, subs = make_item_world(rng, nitem=36, nrep=6, nsub=6, n1d=1, k=k, shared_runs=True, run_offsets=False)
    rows = write_item_world(d, names, subs, rng)
    write_item_table(d / "table.txt", rows, with_cat=True)
    fold_opt = ["-fold_column", "Category"] if folds_by == "category" else []
    run([prog, "-prefix", d / "out", "-axis", "item", "-nfeatures", k, "-algorithm", algo, "-niter", niter, "-tol", 0,
         "-dataTableFile", d / "table.txt", "-cv", mode, "-save_r2", "-progress", "off", "-overwrite"] + fold_opt)
    sub_order, keys, keyed = expected_keys(rows)
    by_label = {s["label"]: s for s in subs}
    pre = [item_preprocess(by_label[lab], keyed[lab], keys)[0] for lab in sub_order]
    n, T = len(pre), len(keys)
    col_of = "cat" if folds_by == "category" else "run"
    level_of_key = [next(r for kk, r in keyed[sub_order[0]] if kk == key)[col_of] for key in keys]
    levels = list(dict.fromkeys(level_of_key))
    folds = [np.array([j for j in range(T) if level_of_key[j] == lv]) for lv in levels]
    head, body = read_scores(d / "out.scores_folds.txt")
    check("scores name the metric 'ident' on the item axis", float(head[5] != "ident"), 0.0, "{:.0f}")
    check("folds are the levels, in order of first appearance", float([r[2] for r in body[:len(levels)]] != levels), 0.0, "{:.0f}")
    got = {(r[0], r[2]): [num(r[4]), num(r[5]), num(r[6]), num(r[7])] for r in body}
    want, vox_ref = {}, [np.zeros(z.shape[0]) for z in pre]
    if mode == "run":
        for fi_, cols in enumerate(folds):
            tr = [np.delete(z, cols, axis=1) for z in pre]
            te = [z[:, cols] for z in pre]
            r = reference_fit(cli, work, f"{name}_f{fi_}", tr, k, niter, algo)
            P = [r["w"][i].T @ te[i] for i in range(n)]
            for i in range(n):
                oth = [j for j in range(n) if j != i]
                ref_u = sum(P[j] for j in oth) / len(oth)
                if algo == "prob":
                    prec = np.linalg.inv(r["sigma"]) + (1.0 / r["rho2"][oth]).sum() * np.eye(k)
                    ref_r = np.linalg.solve(prec, sum(P[j] / r["rho2"][j] for j in oth))
                else:
                    ref_r = ref_u
                acc, ch = np_item_identify(P[i], ref_u)
                sse = ((te[i] - r["w"][i] @ ref_r) ** 2).sum(1)
                sst = ((te[i] - te[i].mean(1, keepdims=True)) ** 2).sum(1)
                want[(sub_order[i], levels[fi_])] = [np_isc(P[i], ref_u), acc, ch, 1.0 - sse.sum() / sst.sum()]
        worst_proj = 0.0
        for i in range(n):
            if by_label[sub_order[i]]["is1D"]:
                continue
            pj = load_afni(d / f"out.proj.{sub_order[i]}").reshape(k, T)
            for fi_, cols in enumerate(folds):
                r = run_ours(cli, work, f"{name}_f{fi_}", algo, [np.delete(z, cols, axis=1).astype(np.float32) for z in pre],
                             k, niter, False, extra=["--canon"])
                worst_proj = max(worst_proj, procrustes_residual(pj[:, cols], r["w"][i].T @ pre[i][:, cols]))
        check("held-out item projections equal the fold models' up to a rotation of the gauge", worst_proj, 2e-3)
    else:
        full = reference_fit(cli, work, name + "_full", pre, k, niter, algo)
        for i in range(n):
            oth = [j for j in range(n) if j != i]
            dd = work / f"{name}_drop{i}"
            dd.mkdir(parents=True, exist_ok=True)
            extra = []
            if algo == "prob":
                full["rho2"][oth].astype(np.float64).tofile(dd / "rho2.f64")
                full["sigma"].astype(np.float32).tofile(dd / "sigma.f32")
                extra = ["--init-rho2", str(dd / "rho2.f64"), "--init-sigma", str(dd / "sigma.f32")]
            r = run_ours(cli, work, f"{name}_drop{i}", algo, [pre[j].astype(np.float32) for j in oth], k, 5, False,
                         init_w=[full["w"][j] for j in oth], extra=extra)
            for fi_, cols in enumerate(folds):
                tr, te = np.delete(pre[i], cols, axis=1), pre[i][:, cols]
                W = np_polar(tr @ np.delete(r["s"], cols, axis=1).T)
                Pm, s_te = W.T @ te, r["s"][:, cols]
                acc, ch = np_item_identify(Pm, s_te)
                sse = ((te - W @ s_te) ** 2).sum(1)
                sst = ((te - te.mean(1, keepdims=True)) ** 2).sum(1)
                want[(sub_order[i], levels[fi_])] = [np_isc(Pm, s_te), acc, ch, 1.0 - sse.sum() / sst.sum()]
    errs = np.array([[abs(got[k_][q] - want[k_][q]) for q in range(4)] for k_ in want])
    check("isc per subject and fold vs the reference", float(errs[:, 0].max()), 5e-4)
    check("item identification accuracy vs the reference", float(errs[:, 1].max()), 1e-6)
    check("... and its chance level (1 over the fold's items)", float(errs[:, 2].max()), 1e-6)
    check("held-out R^2 vs the reference", float(errs[:, 3].max()), 5e-4)
    acc_mean = float(np.mean([want[k_][1] for k_ in want]))
    ch_mean = float(np.mean([want[k_][2] for k_ in want]))
    check("identification is well above chance on this signal (sanity)", float(not acc_mean > 3 * ch_mean), 0.0, "{:.0f}")
    if folds_by == "category":                          # Category varies over a subject's items: not carried along
        head_s, _ = read_scores(d / "out.scores.txt")
        check("an item-level column is not carried into scores.txt", float("Category" in head_s), 0.0, "{:.0f}")
        bad = [dict(r, cat="tool") if (r["subj"] == sub_order[2] and r["item"] == "it01") else r for r in rows]
        write_item_table(d / "badfold.txt", bad, with_cat=True)
        item_fails(prog, d, ["-axis", "item", "-cv", "run", "-fold_column", "Category"], "a fold must hold out the same items",
                   "E: a fold column that differs between subjects for an item", d / "badfold.txt")


def case_time_zscore_all(prog, cli, work, rng):
    """-zscore all on the time axis: each feature z-scored over the whole timeline, run boundaries ignored."""
    print("\n[time axis, -zscore all]", flush=True)
    d = work / "time_zall"
    _, subs, common = make_cv_world(rng, RUNS)
    rows = write_world(d, subs)
    write_table(d / "table.txt", rows, False)
    write_nifti(d / "mask.nii.gz", common.astype(np.float32)[..., None])
    run([prog, "-prefix", d / "out", "-nfeatures", K, "-niter", 40, "-tol", 0, "-zscore", "all", "-mask", d / "mask.nii.gz",
         "-dataTableFile", d / "table.txt", "-progress", "off", "-overwrite"])
    pre = []
    for sub in subs:
        x = sub["x"].astype(np.float32).astype(np.float64)
        x = x[x.std(1) > 1e-7 * np.abs(x.mean(1))]       # only a feature flat over everything is dropped
        pre.append((x - x.mean(1, keepdims=True)) / x.std(1, keepdims=True))
    ref = reference_fit(cli, work, "time_zall_ref", pre, K, 40, "prob")
    check("-zscore all (time axis) vs the driver on data z-scored over the whole timeline",
          relfro(np.loadtxt(d / "out.S.1D").T, ref["s"]), 2e-4)
    hdr, _, _, _ = parse_srm_file(d / "out.srm")
    check("the model records zscore all", float(hdr["zscore"][0] != "all"), 0.0, "{:.0f}")


def write_xmat(path, nfull, starts, good=None):
    """A design matrix in afni_proc.py's X.xmat.1D form (commented NIML), with one column of ones: only the
    header attributes matter here."""
    good = list(range(nfull)) if good is None else good
    rows = "\n".join("1" for _ in good)
    Path(path).write_text("# <matrix\n#  ni_type = \"1*double\"\n"
                          f"#  ni_dimen = \"{len(good)}\"\n#  RowTR = \"2\"\n"
                          f"#  GoodList = \"{','.join(map(str, good))}\"\n#  NRowFull = \"{nfull}\"\n"
                          f"#  RunStart = \"{','.join(map(str, starts))}\"\n#  ColumnLabels = \"Run#1Pol#0\"\n# >\n"
                          f"{rows}\n# </matrix>\n")


def case_concat(prog, cli, work, rng):
    """One input per subject holding every run (afni_proc.py's errts): -concat, or an Xmat column, gives the runs."""
    print("\n[runs inside one input: -concat and Xmat]", flush=True)
    d = work / "concat"
    _, subs, common = make_cv_world(rng, RUNS)
    rows = write_world(d, subs)
    write_table(d / "runs.txt", rows, False)
    write_nifti(d / "mask.nii.gz", common.astype(np.float32)[..., None])
    starts = [0, RUNS[0], RUNS[0] + RUNS[1]]
    T = sum(RUNS)
    one, xm = ["Subj Group InputFile"], ["Subj Group Xmat InputFile"]
    for sub in subs:
        imgs = [nib.load(str(d / f"{sub['label']}_run{r + 1}.nii.gz")) for r in range(len(RUNS))]
        f = d / f"errts.{sub['label']}.nii.gz"
        nib.save(nib.Nifti1Image(np.concatenate([np.asarray(i.dataobj) for i in imgs], axis=3), imgs[0].affine), str(f))
        write_xmat(d / f"X.{sub['label']}.xmat.1D", T, starts)
        one.append(f"{sub['label']} {sub['group']} {f}")
        xm.append(f"{sub['label']} {sub['group']} {d / ('X.' + sub['label'] + '.xmat.1D')} {f}")
    (d / "one.txt").write_text("\n".join(one) + "\n")
    (d / "xmat.txt").write_text("\n".join(xm) + "\n")
    base = [prog, "-nfeatures", K, "-niter", 40, "-tol", 0, "-mask", d / "mask.nii.gz", "-progress", "off", "-overwrite",
            "-cv", "run", "-match_window", 5]
    run(base + ["-prefix", d / "runs", "-dataTableFile", d / "runs.txt"])
    p = run(base + ["-prefix", d / "cat", "-dataTableFile", d / "one.txt", "-concat", "1D: " + " ".join(map(str, starts))])
    log_cat = p.stdout + p.stderr
    p = run(base + ["-prefix", d / "xm", "-dataTableFile", d / "xmat.txt"])
    log_xm = p.stdout + p.stderr
    S0 = np.loadtxt(d / "runs.S.1D")
    check("-concat: S identical to separate run files", float(np.abs(np.loadtxt(d / "cat.S.1D") - S0).max()), 0.0, "{:.2e}")
    check("Xmat column: S identical to separate run files", float(np.abs(np.loadtxt(d / "xm.S.1D") - S0).max()), 0.0, "{:.2e}")
    def scores(f):                                      # every column but the fold's name (r1 from a Run column, 1 here)
        return [r[:2] + r[3:] for r in read_scores(d / f)[1]]
    f0 = scores("runs.scores_folds.txt")
    check("-concat: the -cv run scores are the separate-run ones", float(scores("cat.scores_folds.txt") != f0), 0.0, "{:.0f}")
    check("Xmat: the -cv run scores are the separate-run ones", float(scores("xm.scores_folds.txt") != f0), 0.0, "{:.0f}")
    check("the log reports the runs found in the inputs", float("each input holds 3 runs" not in log_cat or "RunStart" not in log_xm), 0.0, "{:.0f}")
    p = run([prog, "-nfeatures", K, "-niter", 5, "-mask", d / "mask.nii.gz", "-progress", "off", "-overwrite",
             "-prefix", d / "long", "-dataTableFile", d / "one.txt"])
    check("one long input without -concat: the log suggests -concat or Xmat", float("give -concat or an Xmat column" not in p.stdout + p.stderr), 0.0, "{:.0f}")

    # a censored GoodList: reported, not yet applied
    lab0 = subs[0]["label"]
    write_xmat(d / "X.cens.xmat.1D", T, starts, good=[t for t in range(T) if t not in (5, 6, 70)])
    (d / "xmat_cens.txt").write_text((d / "xmat.txt").read_text().replace(f"X.{lab0}.xmat.1D", "X.cens.xmat.1D"))
    p = run(base + ["-prefix", d / "xc", "-dataTableFile", d / "xmat_cens.txt", "-cv", "none"])
    check("a censored GoodList is reported as not applied (yet)", float("not applied in this build" not in p.stdout + p.stderr), 0.0, "{:.0f}")

    def fails(args, needle, label):
        pr = run([prog, "-progress", "off", "-prefix", d / "bad", "-overwrite", "-nfeatures", K, "-mask", d / "mask.nii.gz"] + args,
                 expect_ok=False)
        out = pr.stdout + pr.stderr
        check(label, float(pr.returncode == 0 or needle not in out), 0.0, "{:.0f}")
        if pr.returncode == 0 or needle not in out:
            print(out)
    fails(["-dataTableFile", d / "one.txt", "-concat", "1D: 5 60 110"], "must start at 0", "E: -concat not starting at 0")
    fails(["-dataTableFile", d / "one.txt", "-concat", "1D: 0 110 60"], "must increase", "E: -concat out of order")
    fails(["-dataTableFile", d / "one.txt", "-concat", f"1D: 0 60 {T}"], "the last run starts at", "E: -concat past the end")
    fails(["-dataTableFile", d / "runs.txt", "-concat", "1D: 0 60 110"], "Use one or the other", "E: -concat with a Run column")
    fails(["-dataTableFile", d / "xmat.txt", "-concat", "1D: 0 60 110"], "not both", "E: -concat and an Xmat column")
    fails(["-axis", "item", "-dataTableFile", d / "one.txt", "-concat", "1D: 0 60 110"], "it is for -axis time", "E: -concat with -axis item")
    write_xmat(d / "X.short.xmat.1D", T - 1, starts)
    (d / "xmat_short.txt").write_text((d / "xmat.txt").read_text().replace(f"X.{lab0}.xmat.1D", "X.short.xmat.1D"))
    fails(["-dataTableFile", d / "xmat_short.txt"], "(NRowFull), but the input has", "E: an Xmat whose NRowFull is not the input's length")
    write_xmat(d / "X.other.xmat.1D", T, [0, 90])
    (d / "xmat_other.txt").write_text((d / "xmat.txt").read_text().replace(f"X.{subs[1]['label']}.xmat.1D", "X.other.xmat.1D"))
    fails(["-dataTableFile", d / "xmat_other.txt"], "every subject needs the same runs", "E: Xmat runs that differ between subjects")


def main_items(args, cli, rng):
    case_time_zscore_all(args.prog, cli, args.workdir, rng)
    for algo in ("prob", "det"):
        case_items(args.prog, cli, args.workdir, rng, algo)
    case_items(args.prog, cli, args.workdir, rng, "prob", zmode="run")
    case_items(args.prog, cli, args.workdir, rng, "prob", zmode="none")
    case_items(args.prog, cli, args.workdir, rng, "prob", use_runs=False)
    for mode in ("run", "loso"):
        for algo in ("prob", "det"):
            case_items_cv(args.prog, cli, args.workdir, rng, mode, algo, "run")
    case_items_cv(args.prog, cli, args.workdir, rng, "run", "prob", "category")


def main_cv(args, cli, rng):
    for mode in ("run", "loso"):
        for algo in ("prob", "det"):
            case_cv(args.prog, cli, args.workdir, rng, mode, algo, [60, 50, 70])
    case_cv(args.prog, cli, args.workdir, rng, "run", "prob", [140], label="_blocks")
    case_nfeat_scan(args.prog, cli, args.workdir, rng)


def main_stage3(args, cli, rng):
    for algo in ("prob", "det"):
        case_model_in(args.prog, cli, args.workdir, rng, algo)
    case_apply_W(args.prog, cli, args.workdir, rng, "run")
    case_apply_W(args.prog, cli, args.workdir, rng, "none")
    case_1D_apply_W(args.prog, cli, args.workdir, rng)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prog", type=Path, default=REPO / "build" / "targets_built" / "3dSRM")
    ap.add_argument("--workdir", type=Path, default=Path("check_3dsrm"))
    ap.add_argument("--srm-cli", help="a built srm_cli; compiled from the repo if omitted")
    ap.add_argument("--only-stage3", action="store_true", help="only the -model_in, -update_S and -apply_W cases")
    ap.add_argument("--only-items", action="store_true", help="only the -axis item cases")
    ap.add_argument("--only-concat", action="store_true", help="only the -concat / Xmat cases")
    ap.add_argument("--require-deps", action="store_true", help="fail, not skip, when numpy or nibabel is missing")
    args = ap.parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    cli = build_cli(args.workdir, args.srm_cli)
    rng = np.random.default_rng(12)
    if args.only_concat:
        case_concat(args.prog, cli, args.workdir, rng)
    elif args.only_items:
        main_items(args, cli, rng)
    elif args.only_stage3:
        main_stage3(args, cli, rng)
    else:
        case_basic(args.prog, cli, args.workdir, rng, with_mask_col=False, zcensor=False, algo="prob")
        case_basic(args.prog, cli, args.workdir, rng, with_mask_col=True, zcensor=True, algo="det")
        case_1D_and_errors(args.prog, cli, args.workdir, rng)
        main_items(args, cli, rng)
        case_concat(args.prog, cli, args.workdir, rng)
        main_cv(args, cli, rng)
        main_stage3(args, cli, rng)
    print()
    if failures:
        print(f"FAILED ({len(failures)}): " + "; ".join(failures))
        sys.exit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
