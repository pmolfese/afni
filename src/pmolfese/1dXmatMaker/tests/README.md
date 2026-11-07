<!--
Developed at the NIH

Part of AFNI. See the repository LICENSE.txt for licensing terms and
disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
Third-party and adapted portions retain their existing attribution and license
terms; see THIRD_PARTY.md.
-->

# Testing and packaging 1dXmatMaker

Run from the package root. Rust tests need only the pinned compiler/linker:

```sh
cargo test --frozen
cargo build --release --bins --examples --frozen
```

The scientific suite needs Python 3.10+ with NumPy/SciPy and AFNI consumers:

```sh
python3 tests/run.py suite --maker target/release/1dXmatMaker \
  --afni-bin /absolute/path/to/afni/bin --probes target/release/examples
```

Reports default to `target/validation`; choose `--tag NAME` or `--report-dir DIR`
for another run. `tests/run.py --help` lists focused commands. Add `--help` after
any command for its options. Cargo tests are a separate gate, also run by Make's
`check` and CMake's ctest target.

`fixtures/` holds independent reference values and capture provenance; normal
tests never change them. The layout groups data by reference implementation:

- `afni/compatibility/`: shared reference cases and pinned AFNI executable provenance.
- `afni/baselines/`: polynomial baselines, sampled signals and nuisance.
- `afni/events/`: event timing and event-to-column identities.
- `afni/basis/`: response waveforms and expression functions.
- `afni/semantics/`: modulation, censoring, contrasts, selectors and automatic polort.
- `glmsingle/`: upstream HRF preparation code and expected sampled waveforms.

Each AFNI group has a case catalog and a `reference/` tree with hashes of the
captured files. Baseline/event catalogs use the historical `m0_case` key to reuse
shared compatibility cases without copying their matrices. The merged semantics
manifest retains original capture hashes and case membership; matrix and input
bytes are unchanged. Group names describe behavior rather than development stages.

`tests/validation/source_inventory.json` is the pinned
AFNI interface inventory used by basis coverage checks. Check these without
regenerating expectations:

```sh
python3 tests/run.py reference audit
# Requires the enclosing AFNI source checkout:
python3 tests/run.py inventory --check
```

Reference capture, `capture-hrfs`, and inventory regeneration are explicit
maintenance operations. Review upstream changes before replacing expectations.
Rust integration tests are `tests/*.rs`; API probes live in `tests/probes/`;
Python checks and helpers live in `tests/validation/`. Preserve independent
oracles and documented numerical tolerances when updating them.

## Release archives

Packaging and archive verification require Python 3.12+. Supply a tested native
binary and its actual Rust target triple; packaging does not cross-compile or
qualify a platform. Use a new output directory; existing archives are not replaced.

```sh
python3 tests/package_release.py --binary target/release/1dXmatMaker \
  --target aarch64-apple-darwin --output-dir target/release-packages
python3 tests/run.py package --release-dir target/release-packages \
  --report target/release-packages/verification.json
```

The source archive includes tests, fixtures and compile-time data. The binary
archive includes the executable, README, examples and notices;
it excludes tests, fixtures and the already-embedded waveform table. Sorted
entries and normalized timestamps make archives reproducible for identical input
bytes. Verification builds/tests the extracted source with an empty Cargo cache
and compares example outputs against the supplied binary.

Keep generated CMake files under `target/cmake`. All of `target/` can be deleted;
it is excluded from source archives along with Python bytecode caches.

The workflow gate exercises `1dXmatTool -beta_manifest` through actual MVPA,
classic RSA and runwise RSA fits. The separate AFNI test
`tests/scripts/check_1dXmatTool_betas.py` (relative to the AFNI repository root)
checks fixed-design Picard mappings, annotations, quoting and rejection cases;
its `--help` lists the explicit executable paths required.
