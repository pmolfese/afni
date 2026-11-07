<!--
Developed at the NIH

Part of AFNI. See the repository LICENSE.txt for licensing terms and
disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
Third-party and adapted portions retain their existing attribution and license
terms; see THIRD_PARTY.md.
-->

# Scientific source attribution and redistribution notes

This package is developed inside AFNI and remains subject to the repository's
[AFNI license](licenses/AFNI.txt). No external Rust crate, C/Fortran object,
AFNI shared library, or native build script is linked into this package.

The response functions, endpoint/amplitude conventions, spline/polynomial
expressions, WAV transition, and duration-table construction in `basis.rs` and
`basis_extended.rs` adapt AFNI `3dDeconvolve.c`. That file credits B. Douglas Ward
and the Medical College of Wisconsin (copyright 1998–2003). AFNI's major MCW
portions are released under Creative Commons Attribution 4.0 (CC BY 4.0):
<https://creativecommons.org/licenses/by/4.0/>. Preserve this attribution and
the upstream license when redistributing this adaptation.

Gamma peak/FWHM conversion and its bounded search adapt `cs_gamfit.c` and
`powell_int.c` within AFNI. Expression language/evaluation conventions and
specific helpers adapt `parser.c`, `parser.inc`, `parser_int.c`, `mri_stats.c`
and `testprime.c`; the parser/statistics files carry the MCW 1994–2000 notice.
All source references are to revision
2664fdd5de7425f658c6d901765fbb4c291693fa. The numerical routines were rewritten
in Rust with owned local state, checked parsing and explicit error handling.
The BLOCK integral was algebraically simplified; probability helpers use
independently expressed mathematical series/continued fractions and bracketing.
The random-stream policy differs from AFNI and is documented in [README.md](README.md#basis-functions-and-expressions).

NIfTI intent codes and distribution/parameter conventions follow
`src/nifti/nifticdf/nifticdf.c`. Its header credits Brown and Lovato's public-domain
cdflib, K. Krishnamoorthy's unrestricted-use code, and RW Cox's public-domain
portions. No cdflib/NIfTI implementation is bundled or linked here. The header
also disclaims warranty or institutional liability; see that source and
`src/nifti/LICENSE` for its full notice. These provenance notes do not imply
institutional endorsement of the Rust implementation.

Frozen small synthetic waveform/matrix fixtures were generated locally by the
identified AFNI executable. Their capture commands, inputs, executable hash and
output hashes accompany them. No participant data or external research dataset
is redistributed. SciPy is used only by an optional validation script through
the developer's installation; its source/binaries are not vendored or linked.

## GLMsingle sampled library

The embedded `data/glmsingle/getcanonicalhrflibrary.tsv` and unmodified reference
Python files at `glmsingle/upstream/` inside `fixtures.tgz` come from
[cvnlab/GLMsingle](https://github.com/cvnlab/GLMsingle/tree/082851f9db87ff30c49370f27365ac1112f18ff7),
commit **082851f9db87ff30c49370f27365ac1112f18ff7**. Copyright (c) 2021,
Kendrick Kay; licensed under BSD-3-Clause. Preserve the complete notice and
conditions in [data/glmsingle/LICENSE](data/glmsingle/LICENSE). The exact source
hashes, sample spacing, original stimulus duration and shape are recorded in
[data/glmsingle/provenance.json](data/glmsingle/provenance.json).

`sampled.rs` independently implements the pinned duration summation, PCHIP
resampling and float32 normalization sequence. The PCHIP slopes follow the
Fritsch–Butland mathematical convention used by the upstream SciPy call; no
SciPy implementation is bundled or linked. Frozen fixtures were captured by
running the unmodified pinned Python library with NumPy/SciPy. Its external
validation dependencies do not become runtime dependencies of the Rust program.

Source and binary distributions include AFNI and NIfTI notices in `licenses/`;
binary installations also include the GLMsingle license as
`share/1dXmatMaker/licenses/GLMsingle.txt`. These notices preserve scientific and
copyright attribution, and do not imply endorsement or complete GLMsingle fitter
compatibility. See [README.md](README.md#sampled-hrf-candidates) for the precise scope.
