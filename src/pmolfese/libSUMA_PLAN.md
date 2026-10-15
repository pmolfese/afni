# Plan: split libSUMA into headless algorithms and GUI support

## Executive summary

This split is possible, but it is not currently a one-line CMake source-list change. The main obstacle is API coupling, not the linker:

- `SUMA_suma.h` unconditionally includes `SUMA_X_objs.h`, which includes X11, Xt, Motif, OpenGL, GLX, and X-dependent `mrix` headers.
- `SUMA_define.h` mixes computational data structures with widget, X11, and OpenGL types. A scan of that header finds about 234 uses of types such as `Widget`, `Display`, `XEvent`, `GLfloat`, and `GLXContext`.
- Several large source files mix algorithms/object lifetime code with display code. In particular, `SUMA_Color.c`, `SUMA_CreateDO.c`, `SUMA_DOmanip.c`, `SUMA_GeomComp.c`, `SUMA_Load_Surface_Object.c`, `SUMA_Surface_IO.c`, and `SUMA_niml.c` cannot be assigned cleanly to one library without moving functions or introducing a narrow interface.
- The current CMake `SUMA` target contains 44 C translation units and conditionally adds the bundled GLw implementation. It publicly exposes GLUT/GLw and privately links `AFNI::mrix`, GLU, Xmu, and other dependencies. Consequently, command-line programs such as `MapIcosahedron`, `SurfClust`, `SurfSmooth`, and `ConvertSurface` acquire X/Motif/GL dependencies even when they do not open a window.

## Fixed decisions and scope

The following are requirements, not open questions:

1. Compatibility means **source and behavior compatibility**, not binary ABI compatibility. Existing `3d*` and `Surf*` programs must retain their command-line interfaces, accepted inputs, output formats/names, exit behavior, and numerical results. Their sources may be adjusted to include/link the new core API.
2. `libSUMA` remains the full compatibility library containing both the headless implementation and the GUI implementation. Existing full-SUMA consumers continue to use `-lSUMA` or `AFNI::SUMA`.
3. `libSUMA_core` is the new headless library. It must contain no X11, Xt, Xmu, Motif, OpenGL, GLU, GLX, GLUT, GLw, `mrix`, or `coxplot` dependency.
4. Both CMake and the legacy Makefile build remain supported and must receive the split in the same change series.
5. All supported `3d*` and `Surf*` programs must continue to build and behave the same. Programs that only need surface algorithms should ultimately link the core library; graphical programs may continue to link full `libSUMA`.
6. “Pure algorithms” means free of graphical dependencies and UI side effects. AFNI MRI/NIML, NIfTI, GSL, GTS, zlib, OpenMP, and other non-graphical data/algorithm libraries remain allowed when needed.

The desired installed artifacts and consumer paths are:

```text
                              shared core object set
                                      │
                         ┌────────────┴────────────┐
                         │                         │
headless programs ─> libSUMA_core            libSUMA <─ GUI/full consumers
                     AFNI::SUMA_core          AFNI::SUMA
                         │                         │
                         └─> non-graphical deps   ├─> same core implementation
                                                   └─> X/Motif/GL/mrix/etc.

optional CMake alias: AFNI::SUMA_gui ──────────────> AFNI::SUMA
```

Core code must never call GUI code. Full `libSUMA` is assembled from the core object set plus GUI sources so it continues to define the complete historical symbol set itself. This is more compatible with direct legacy `-lSUMA` links than making `libSUMA` depend transitively on `libSUMA_core`. No executable should link both libraries.

## Proposed names and compatibility policy

Use these CMake targets/artifacts:

- `SUMA_core`, exported as `AFNI::SUMA_core`, producing `libSUMA_core`.
- `SUMA`, exported as `AFNI::SUMA`, continuing to produce `libSUMA` and containing the core plus X/Motif/OpenGL implementation.
- Optionally add `AFNI::SUMA_gui` as an alias of `SUMA` to make intent clear in new CMake code. Do not produce a third installed `libSUMA_gui` artifact.

Compile the core source list once as position-independent objects, then use those objects in both installed libraries:

```cmake
add_library(SUMA_core_objects OBJECT ${SUMA_CORE_SOURCES})
set_property(TARGET SUMA_core_objects PROPERTY POSITION_INDEPENDENT_CODE ON)

add_afni_library(SUMA_core $<TARGET_OBJECTS:SUMA_core_objects>)
add_afni_library(SUMA
  $<TARGET_OBJECTS:SUMA_core_objects>
  ${SUMA_GUI_SOURCES})
add_library(AFNI::SUMA_gui ALIAS SUMA)
```

Attach non-graphical usage requirements to both installed targets and graphical requirements only to `SUMA`. A small helper function can avoid duplicating the core link-dependency declaration. Do **not** link `SUMA` to `SUMA_core` when it already embeds the core objects; that risks duplicate definitions and changes legacy shared-library resolution behavior.

`SUMA_core` must remain a real, independently linkable installed artifact. The internal object target is only a way to keep the two libraries' implementation identical.

## Dependency rules

### Allowed in `SUMA_core`

The exact list should be minimized while the split is implemented, but these are reasonable non-GUI dependencies:

- C standard library and system math library
- `AFNI::mri` (not `AFNI::mrix`)
- AFNI's non-GUI NIML support
- NIfTI libraries
- `AFNI::track_tools`, `AFNI::segtools`, `AFNI::3DEdge`, GSL, and GTS only where an extracted algorithm actually needs them
- zlib and OpenMP where already required

Every dependency should be `PRIVATE` unless a type or required definition genuinely appears in the core public API.

### Forbidden in `SUMA_core`

- `AFNI::mrix` and `AFNI::coxplot`
- X11, Xt, Xmu, Motif/Xm, GLw
- OpenGL, GLU, GLX, and GLUT
- `imseq.h`, `display.h`, `xim.h`, `xutil.h`, and other X-facing headers
- widget handles or GL scalar/object types in public core structures and function signatures
- event-pump, popup, redisplay, snapshot, or window-management calls from core code

Enforce this mechanically in CI; do not rely on target naming or developer convention.

## Source ownership map

This is a starting map, not a promise that each whole file can move unchanged. Confirm each item with compilation and undefined-symbol checks as the work proceeds.

### Strong initial candidates for `SUMA_core`

- `MarchingCubes/MarchingCubes.c`
- `MarchingCubes/mc_ziad.c`
- `PLY/obj_io.c`
- `PLY/plyfile.c`
- `SUMA_DataSets.c`
- `SUMA_ExpEval.c` after moving `SUMA_DispExpr_To_EvalExpr` to core
- `SUMA_HomerFunc.c`
- `SUMA_IsoSurfaceFunc.c`
- `SUMA_IV_FaceSetsextract.c`
- `SUMA_IV_XYZextract.c`
- `SUMA_LocalStat.c`
- `SUMA_MiscFunc.c`
- `SUMA_SegFunc.c`
- `SUMA_SegOpts.c`
- `SUMA_spharm.c`
- `SUMA_SphericalMapping.c`
- `SUMA_SurfaceToSurface.c`
- `SUMA_SurfNorm.c`
- `SUMA_trackball.c`
- `SUMA_VolData.c`

These files still need their `SUMA_suma.h` includes replaced by narrow core headers. Some depend on allocation, surface-object, or colormap functions currently living in mixed files; those supporting functions must move to core first.

### Strong candidates for `SUMA_gui`

- `SUMA_display.c`
- `SUMA_display_b.c`
- `SUMA_input.c`
- `SUMA_SVmanip.c`
- GUI/event portions of `SUMA_Engine.c`
- `SUMA_plot.c`
- `SUMA_volume_render.c`
- `SUMA_xColBar.c`
- `GLw_local/GL/GLwDrawA.c`
- the SUMA GUI entry point, `SUMA_suma.c`

The current object files confirm direct X/Motif/GL references in these areas. `SUMA_Color.c` and `SUMA_CreateDO.c` also contain many such references, but they own important non-GUI functions and therefore belong in the mixed category below.

### Mixed files that should be split by responsibility

| Current file | Core responsibility to extract | GUI responsibility to retain/extract |
|---|---|---|
| `SUMA_Color.c` | colormap parsing, lookup, interpolation, scaling, color-vector allocation | X color allocation, widget synchronization, OpenGL material/drawing helpers |
| `SUMA_CreateDO.c` | allocation/free/accessors for surfaces and other displayable data objects; dimensions and metadata | OpenGL drawing, texture creation, viewer/controller setup |
| `SUMA_DOmanip.c` | data-object registry, accessors, relationships, nonvisual transforms | window raising, controller updates, redisplay and UI actions |
| `SUMA_GeomComp.c` | geometry and mesh algorithms | calls that trigger redisplay or manipulate viewer state |
| `SUMA_Load_Surface_Object.c` | file parsing and construction of surface/data objects | attaching overlays/controllers and updating viewers |
| `SUMA_Surface_IO.c` | surface/ROI serialization and format conversion | dialogs, forced user questions, controller initialization, focus changes |
| `SUMA_BrainWrap.c` | numerical brain-wrap algorithm | progress window/event pumping and fallback resources |
| `SUMA_niml.c` | protocol serialization/deserialization and data conversion | command dispatch into viewer, ROI UI, redisplay, work-proc scheduling |
| `SUMA_ParseCommands.c` | command parsing and validation | popup messages, engine dispatch, display-only color helpers |
| `SUMA_help.c` | static help text and source/version information | opening help windows and snapshot-related UI helpers |
| `SUMA_driver.c` | headless driver parsing/data operations, if reusable | Xt application setup, event processing, pause-for-user behavior |
| `SUMA_clippingPlanes.c` | clipping-plane math and serialization | viewer focus, display-object registration, title/redisplay updates |
| `SUMA_dot.c` | numerical dot/overlay computation | overlay installation, controller refresh, redisplay |
| `SUMA_global.c` | only genuinely headless global state, preferably reduced | viewer/controller globals and X resources |

Use responsibility-oriented filenames when extracting code, for example `SUMA_Color_core.c`, `SUMA_Object_core.c`, `SUMA_Surface_IO_core.c`, and `SUMA_niml_core.c`. Avoid creating parallel files that still depend on the umbrella GUI header.

## Header/API refactor

This is the critical path and should precede the final target split.

### 1. Introduce a headless umbrella header

Add `SUMA_core.h` containing only the public headless API. It may aggregate narrowly scoped headers, but it must compile in a translation unit without X11 or OpenGL development headers installed.

Initially it should include only:

- basic SUMA enums/constants
- surface, mesh, dataset, ROI, colormap-data, and algorithm option structures
- core function declarations
- non-GUI AFNI/NIML types that are intentionally part of the API

Keep `SUMA_suma.h` as the compatibility/full API header. It should eventually become conceptually:

```c
#include "SUMA_core.h"
#include "SUMA_gui.h"
```

### 2. Split `SUMA_define.h`

Move definitions into at least two layers:

- `SUMA_core_types.h`: booleans, enums, meshes, surfaces, datasets, ROI data, colormap data, geometry, algorithm options, and nonvisual object metadata.
- `SUMA_gui_types.h`: widgets, callbacks using Xt types, controllers, viewer windows, GL contexts, render records, X resources, and other UI state.

Where a core structure currently contains a GUI member, use one of these approaches in order of preference:

1. Move the GUI member into a separate GUI-side structure keyed by object/viewer ID.
2. Store an opaque pointer declared as an incomplete type, with ownership entirely in `SUMA_gui`.
3. As a temporary migration bridge only, use a `void *gui_private` field and replace it with a typed opaque handle before declaring the split complete.

Do not substitute OpenGL types with typedefs in the core header merely to hide the includes. For example, core coordinates/colors should use ordinary `float`/`double`; GL conversion belongs at the renderer boundary.

### 3. Split declarations by ownership

Create narrow headers such as:

- `SUMA_surface_core.h`
- `SUMA_dataset_core.h`
- `SUMA_color_core.h`
- `SUMA_geometry_core.h`
- `SUMA_io_core.h`
- `SUMA_niml_core.h`
- `SUMA_gui.h`

Move declarations out of the catch-all `SUMA_prototype.h` over time. A source file in `SUMA_core` should include the narrow headers it uses, not `SUMA_suma.h`.

### 4. Remove core-to-GUI callbacks

There are current core-looking objects that reference GUI-owned functions such as `SUMA_handleRedisplay`, `SUMA_Remixedisplay`, `SUMA_PopUpMessage`, controller initialization, and viewer focus. Refactor these into one of:

- a pure function that returns status/change information and lets the GUI react;
- an explicit optional callback interface passed by the GUI;
- a GUI wrapper that calls the core operation and then performs redisplay/UI updates.

Prefer return values and GUI wrappers over a global callback table. They keep the core deterministic and easier to test.

### 5. Establish ownership boundaries

Core must own and free all core objects. GUI may attach GUI-private state but must not redefine core object lifetime. Put allocation, free, copy, and accessors for `SUMA_SURFACE_OBJECT` and other headless data structures in `SUMA_core`, even if those functions currently reside in `SUMA_CreateDO.c`.

## CMake and build-option changes

### Phase A: add a separately buildable core option

Introduce `COMP_SUMA_CORE`, defaulting on when either core SUMA command-line tools or full SUMA are requested. Keep `COMP_SUMA` as the existing GUI-suite option for compatibility.

Desired option relationship:

```text
COMP_SUMA_CORE: requires core libraries, but not COMP_GUI
COMP_SUMA:      requires COMP_SUMA_CORE and COMP_GUI
```

Then:

- call `add_subdirectory(SUMA)` when `COMP_SUMA_CORE OR COMP_SUMA`;
- discover GSL/NIfTI/core dependencies for `COMP_SUMA_CORE`;
- discover X11, Motif, OpenGL/GLX, GLU, GLUT, GLw, and `mrix` only for `COMP_SUMA`/`COMP_GUI`;
- ensure `src/CMakeLists_x_dependent.txt` is not needed for a core-only configuration;
- move GTS discovery/building to the core side if `suma_gts` proves headless after inspection.

### Phase B: define two targets

The final source lists follow the completed file split, but the target construction should have this shape:

```cmake
add_library(SUMA_core_objects OBJECT ${SUMA_CORE_SOURCES})
set_property(TARGET SUMA_core_objects PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(
  SUMA_core_objects PRIVATE PLY ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(SUMA_core_objects PRIVATE SUMA_CORE_BUILD)

add_afni_library(SUMA_core $<TARGET_OBJECTS:SUMA_core_objects>)
target_link_libraries(SUMA_core PRIVATE AFNI::mri ...non-graphical deps...)

if(COMP_SUMA)
  add_afni_library(
    SUMA
    $<TARGET_OBJECTS:SUMA_core_objects>
    ${SUMA_GUI_SOURCES})
  add_library(AFNI::SUMA_gui ALIAS SUMA)
  target_link_libraries(
    SUMA
    PRIVATE ...non-graphical deps...
            AFNI::mrix Motif::Motif X11::Xt X11::X11 X11::Xmu
            ${_SUMA_GLU_TARGETS} ...)
endif()
```

The full target intentionally embeds the core objects rather than linking `SUMA_core`. This ensures `libSUMA` itself retains the historical functions moved into the core source set, including for legacy Make consumers that specify only `-lSUMA`. Add a test that rejects any executable link line containing both `SUMA` and `SUMA_core`.

`SUMA_COMPILED` is currently `PUBLIC`. Audit what it changes and replace it with target-specific definitions (`SUMA_CORE_BUILD` and `SUMA_GUI_BUILD`) or make it `PRIVATE`. Consumers should not need internal build macros merely to include public headers.

### Phase C: classify executables

Split the large `PROGRAMS` loop into explicit `SUMA_CORE_PROGRAMS` and `SUMA_GUI_PROGRAMS` lists. This classification must be identical in purpose in CMake and Make; if maintaining two literal manifests becomes error-prone, store the lists in a simple included `.mk`/`.cmake` generated from a checked-in neutral manifest.

Likely early core-program candidates include:

- `MapIcosahedron`
- `CreateIcosahedron`
- `ConvertSurface`
- `SurfaceMetrics`
- `SurfClust`
- `SurfClustSim`
- `SurfSmooth`
- `SurfInfo`
- `SurfDsetInfo`
- `3dVol2Surf`
- `3dSurf2Vol`
- `SurfToSurf`
- `SurfQual`
- `SurfFWHM`
- `SurfExtrema`
- `SpharmDeco` and `SpharmReco`

This list must be proven one target at a time. A command-line interface is not automatically headless: inspect its undefined symbols and runtime paths after relinking. GUI utilities (`suma`, `prompt_popup`, `SUMA_paperplane`, `SUMA_glxdino`, `SUMA_pixmap2eps`, and interactive/controller tools) remain on `AFNI::SUMA`.

Update `3dRSA` to link `AFNI::SUMA_core` if its surface support uses only the headless API. That would remove the current special-case requirement to link `AFNI::mrix` solely because `SUMA_suma.h` exposes GUI headers.

### Phase D: implement legacy Makefile parity

The source has supported `SUMA_Makefile_NoDev`, `SUMA_Makefile`, and `Makefile.INCLUDE` rules for `libSUMA.a`/`libSUMA.so`; they must be updated, not retired.

- Replace `SUMA_SRCS` with `SUMA_CORE_SRCS` and `SUMA_GUI_SRCS`; derive `SUMA_CORE_OBJS` and `SUMA_GUI_OBJS` mechanically.
- Add a `LIBSUMA_CORE_OBJ` selector beside `LIBSUMA_OBJ` in `Makefile.INCLUDE`, choosing `libSUMA_core.a` or the platform's shared-library form consistently with the existing `MRI_SHARED` logic.
- Build `libSUMA_core.a`/`.so` from `SUMA_CORE_OBJS`, Marching Cubes, and PLY objects only.
- Build compatibility `libSUMA.a`/`.so` from the same `SUMA_CORE_OBJS` plus `SUMA_GUI_OBJS` and GLw objects as required. It must continue to contain all legacy symbols and must not require consumers to add `-lSUMA_core`.
- Split `SUMA_PROGS`/`SUMA_EXECS` into `SUMA_CORE_PROGS` and `SUMA_GUI_PROGS`. Core program recipes use `LIBSUMA_CORE_OBJ` plus non-graphical libraries; GUI recipes retain `LIBSUMA_OBJ` and graphical libraries.
- Add top-level rules that copy/install both libraries. Keep the installed binary names and current `SUMA_PROGS` union unchanged.
- Update `clean`, `anew`, `all`, `standalones`, `install`, and `packbins` so both library artifacts and every current program are covered.
- Preserve the special `SUMA_xColBar.o` optimization rule on the GUI side.
- Exercise both static and shared Make configurations; archive member ordering and shared link flags must be tested separately.

The CMake and Make source manifests must be compared automatically in CI so a later source cannot be added to one library in one build system and a different library in the other.

### Concrete Make target shape

The exact platform suffix continues to come from `MRI_SHARED`, but the rules should reduce to this model:

```make
SUMA_CORE_OBJS = $(SUMA_CORE_SRCS:.c=.o) MarchingCubes.o mc_ziad.o
SUMA_GUI_OBJS  = $(SUMA_GUI_SRCS:.c=.o) $(GLW_OBJS)

libSUMA_core.a: libply.a $(SUMA_CORE_OBJS)
	$(AR) q $@ $(SUMA_CORE_OBJS)
	$(RANLIB) $@

libSUMA.a: libply.a libglut.a $(SUMA_CORE_OBJS) $(SUMA_GUI_OBJS)
	$(AR) q $@ $(SUMA_CORE_OBJS) $(SUMA_GUI_OBJS)
	$(RANLIB) $@
```

Provide parallel shared-library rules using the existing platform flags. Define separate link bundles:

- `SUMA_CORE_LIBLIST`: `LIBSUMA_CORE_OBJ`, `libply`, MRI/NIML/NIfTI, track/seg/GTS, math, and other non-graphical libraries actually needed.
- `SUMA_GUI_LIBLIST`: `LIBSUMA_OBJ` plus GLUT/GL/GLU/GLw, `mrix`/`coxplot`, X11/Motif, and the non-graphical dependencies.

Then change each program recipe to use exactly one bundle. This is a real Makefile refactor; preserving the Make build does not mean preserving its current variables or recipes.

### File-by-file build checklist

- `cmake/afni_cmake_build_options.cmake`: add `COMP_SUMA_CORE`; make `COMP_SUMA` depend on both `COMP_SUMA_CORE` and `COMP_GUI`.
- `cmake/afni_project_dependencies.cmake`: move non-graphical SUMA discovery under `COMP_SUMA_CORE`; leave X/Motif/GL discovery under the GUI/full conditions.
- `src/CMakeLists.txt`: enter `src/SUMA` for core-only builds; link `3dRSA` and `3dMVPA` to the core when enabled.
- `src/SUMA/CMakeLists.txt`: define the core object set, both installed libraries, two program manifests, install/export behavior, and dependency checks.
- `src/Makefile.INCLUDE`: add `LIBSUMA_CORE_OBJ`, top-level build/copy/install dependencies, and core/full program manifests.
- `src/SUMA/SUMA_Makefile_NoDev`: add core/full source and object lists, static/shared rules, link bundles, program dependencies, clean/install/package handling.
- `src/SUMA/SUMA_Makefile`: pass through or mirror the new variables and platform settings used by `SUMA_Makefile_NoDev`.
- `tests/CMakeLists.txt` and `tests/scripts/`: add the header, dependency, manifest, CLI, numerical, and GUI regression gates.
- `packaging/installation_components.txt` and related packaging logic: place `libSUMA_core` and headless programs outside the graphical runtime dependency closure while retaining full `libSUMA` in the SUMA/GUI component.

## Required program inventory

The broad compatibility rule applies to every AFNI program whose installed name begins with `3d` or `Surf`, even if it does not use SUMA. The split directly changes only the programs that consume SUMA code, so use the following checked-in inventory as the focused migration matrix. Preserve exact executable capitalization, including `3dinfill` and `SurfLocalstat`.

### Required headless `3d*` consumers

These programs must build and run against `libSUMA_core` without graphical libraries:

- `3dBRAIN_VOYAGERtoAFNI`
- `3dCRUISEtoAFNI`
- `3dGenFeatureDist`
- `3dGenPriors`
- `3dHist`
- `3dSeg`
- `3dSkullStrip`
- `3dSurf2Vol`
- `3dSurfMask`
- `3dVol2Surf`
- `3dinfill`
- `3dRSA`
- `3dMVPA`

The Makefile also exposes non-default/special targets `3dSeg_noomp` and `3dfilter`; keep both working against the core even though they are not in the default program list.

### Required headless `Surf*` consumers

These programs must build and run against `libSUMA_core` without graphical libraries:

- `Surf2VolCoord`
- `SurfClust`
- `SurfClustSim`
- `SurfDist`
- `SurfDsetInfo`
- `SurfExtrema`
- `SurfFWHM`
- `SurfInfo`
- `SurfLocalstat` (currently present in the Make build; preserve this spelling)
- `SurfMeasures`
- `SurfMesh`
- `SurfPatch`
- `SurfQual`
- `SurfRetinoMap`
- `SurfSmooth`
- `SurfToSurf`
- `SurfaceMetrics`

`SurfMatch` is a named Make target even though it is not in the current default `SUMA_PROGS` list; keep it buildable against the core/GTS stack. `SUMA_SurfNorm.c`, `SUMA_SurfWarp.c`, `SUMA_SurfaceToSurface.c`, and `SUMA_Surface_IO.c` are implementation sources, not additional installed program names.

### Manifest reconciliation required before code movement

The current CMake and Make program sets are already not identical (`SurfLocalstat` is a visible example). PR 1 must create a table with one row per program and these columns:

| Program | CMake target | Make target | Installed today | Intended library | Smoke command | Numerical/golden test |
|---|---:|---:|---:|---|---|---|
| Example: `SurfClust` | yes | yes | yes | `SUMA_core` | `SurfClust -help` | existing/new fixture |

The checked-in table becomes the authoritative migration checklist. It should also include all remaining `3d*` and `Surf*` programs outside SUMA with `Intended library = unaffected`, so the broad regression promise is measurable rather than implicit.

## Implementation sequence

Keep each step buildable and behavior-preserving.

### PR 1: dependency audit and guardrails

- Add the required program inventory described above, reconciling `src/SUMA/CMakeLists.txt`, `SUMA_PROGS`/`SUMA_EXECS`, named Make targets, installation rules, and packaging rules.
- Add a script or CMake test that fails if the candidate core target links forbidden libraries.
- Add a compile-only public-header test for `SUMA_core.h` in an environment/include path without X/GL headers.
- Record the current exported symbols of `libSUMA` (`nm`) and current runtime dependencies (`ldd`/`otool -L`) for comparison.
- Capture baseline `-help` output/exit status and at least one meaningful input/output fixture for every directly affected `3d*` and `Surf*` program. Mark programs lacking a usable fixture; do not move them until one exists.
- Add a Make/CMake manifest comparison check.

Deliverable: no production behavior change, but objective checks exist.

### PR 2: separate headers and types

- Introduce `SUMA_core_types.h` and `SUMA_core.h`.
- Move pure definitions first; keep compatibility includes/typedefs in `SUMA_suma.h` as needed.
- Replace GL scalar types in core structures with C scalar types.
- Isolate GUI-only structures in `SUMA_gui_types.h`.
- Convert core candidate source files from `SUMA_suma.h` to narrow headers.

Deliverable: a trivial C file including only `SUMA_core.h` compiles without X/Motif/GL include directories.

### PR 3: extract mixed-file core services

- Extract object allocation/lifetime/accessors from `SUMA_CreateDO.c` and `SUMA_DOmanip.c`.
- Extract pure colormap operations from `SUMA_Color.c`.
- Separate surface parsing/construction from viewer attachment in `SUMA_Load_Surface_Object.c` and `SUMA_Surface_IO.c`.
- Make geometry, clipping, NIML, and command parsing return changes/errors instead of directly requesting redisplay or popup dialogs.
- Move progress UI/event pumping out of `SUMA_BrainWrap.c` into a GUI wrapper or optional caller callback.

Deliverable: all strong core candidates form a closed symbol graph with no references into GUI objects.

### PR 4: create both library layouts in CMake and Make

- Add the real `SUMA_core_objects`, `SUMA_core`, and full `SUMA` CMake targets.
- Add `SUMA_CORE_SRCS`/`SUMA_GUI_SRCS`, `libSUMA_core`, and full `libSUMA` rules to the Make build.
- Initially migrate `SurfaceMetrics`, `SurfClust`, and `3dVol2Surf`; together they exercise mesh algorithms, clustering, surface I/O, and volume/surface integration.
- Verify that full `libSUMA` still contains the baseline legacy symbol set and that the `suma` executable still links only full `libSUMA`.
- Compare the three migrated programs against baseline fixtures under both build systems.

Deliverable: the two installed libraries exist under both build systems; the three migrated binaries have no graphical runtime dependency and unchanged behavior; full `libSUMA` still supports the GUI.

### PR 5: migrate every required `3d*` and `Surf*` consumer

- Move the required program inventory incrementally, keeping the CMake and Make rows synchronized.
- Move `suma_gts` to the core dependency stack for `SurfMesh` and `SurfMatch`.
- Migrate `3dRSA` and `3dMVPA` from the full SUMA/mrix pairing to `AFNI::SUMA_core`; retain `USE_SUMA` behavior and remove the special `AFNI::mrix` link when no GUI symbol remains.
- Remove redundant direct dependencies from each executable; let target usage requirements carry only what is necessary.
- After each batch, run baseline CLI/output comparisons and the forbidden-dependency check.

Deliverable: every required `3d*`/`Surf*` consumer builds under CMake and Make, runs in a headless environment, and passes its compatibility fixture.

### PR 6: packaging, compatibility, and cleanup

- Install/export `AFNI::SUMA_core` and its public headers.
- Preserve `AFNI::SUMA` and `libSUMA` permanently as the full compatibility library unless a future, separate deprecation proposal is accepted.
- Update component packaging so core tools do not pull the GUI/SUMA runtime package.
- Verify Make `install`/`packbins` and CMake component installation contain the same required executables and the appropriate libraries.
- Document the new API and rules for placing new code.

Deliverable: downstream users can choose the headless core or the full GUI library without accidental dependency leakage.

## Verification strategy

### Configure/build matrix

At minimum test:

1. Core-only: `COMP_GUI=OFF`, `COMP_SUMA_CORE=ON`, `COMP_SUMA=OFF`.
2. Full GUI/SUMA: existing default GUI configuration.
3. Shared libraries on Linux.
4. Static libraries on Linux, because static archives expose dependency/order mistakes that shared builds can hide.
5. macOS/XQuartz shared build.
6. Bundled versus system GLUT/GLw where those paths remain supported.

### Mechanical dependency checks

For `libSUMA_core` and every migrated headless executable:

- Linux: inspect `readelf -d` or `ldd` and reject X11, Xt, Xm, Xmu, GL, GLU, GLUT, GLX, GLw, `mrix`, and `coxplot`.
- macOS: inspect `otool -L` with the same reject list.
- Inspect undefined symbols with `nm -u`; reject `X*`, `Xt*`, `Xm*`, `Xmu*`, `gl*`, `glu*`, `glut*`, `GLw*`, `ISQ_*`, and GUI-only `MCW_*` symbols. Maintain a small allowlist only if a prefix has both GUI and non-GUI functions.
- Preprocess or compile `SUMA_core.h` without GUI include paths and verify no GUI header enters the include trace.
- Add a source-level check forbidding `SUMA_suma.h` and GUI headers in `SUMA_core` sources.

### Behavior/regression checks

- Run existing `test_SurfClustSim.py` and relevant command tests.
- Add golden/numeric tests for surface loading, mesh neighbors/normals, smoothing, clustering, mapping, colormap scaling, NIML round trips, and surface format conversion.
- For every program in the required inventory, compare the pre-split and post-split `-help` output, exit status, output file names, output headers/metadata, and output data. Use byte-for-byte comparison for deterministic text/binary formats and documented numerical tolerances only where output is already nondeterministic or floating-point-order dependent.
- Run the full existing test selection for all other repository `3d*` and `Surf*` programs, even when their link graph is unaffected.
- Run `tests/scripts/test_guis.py` for the full build to prove that moving core functions did not break `suma` behavior.
- Add allocation/free stress tests or sanitizer runs around surface and overlay lifetime boundaries; this is where splitting `SUMA_CreateDO.c` is most likely to introduce regressions.

### Source/link compatibility checks

- Compare the old `libSUMA` exported-symbol list with the compatibility `libSUMA` after each extraction. All symbols used by in-tree consumers must remain directly defined by full `libSUMA`.
- The exported symbols of `libSUMA_core` are expected to be a subset of the symbols directly exported by full `libSUMA`, because both artifacts contain the same core objects. Reject any unexpected GUI-only symbol in the core, rather than rejecting the intentional overlap.
- Compile every in-tree consumer from source against the revised full `SUMA_suma.h`/`libSUMA` pair and every migrated headless consumer against `SUMA_core.h`/`libSUMA_core`.
- Direct legacy Make links using only `-lSUMA` must continue to resolve without adding `-lSUMA_core`.
- Binary ABI and public structure layout stability are not requirements for this project. Structure changes are allowed when all in-tree callers are updated and behavior tests pass. Still prefer opaque GUI-private state because it creates a durable boundary and reduces future source churn.

## Completion criteria

The split is complete when all of the following are true:

- `SUMA_core.h` compiles without X11, Motif, or OpenGL headers available.
- `libSUMA_core` and every required headless `3d*`/`Surf*` consumer build with `COMP_GUI=OFF`.
- `libSUMA_core` and its headless consumers have no X/Motif/GL/`mrix` runtime dependencies or undefined symbols.
- Full `libSUMA` is built from the core objects plus GUI objects, retains the complete in-tree legacy symbol set, and does not require a second `-lSUMA_core` argument.
- No target links both `libSUMA` and `libSUMA_core`.
- Every required `3d*` and `Surf*` compatibility fixture passes with unchanged CLI and results.
- The full `suma` GUI and GUI test suite still pass.
- Installed CMake exports and package components express the split correctly.
- `AFNI::SUMA`/`libSUMA` remain the documented full compatibility choice; `AFNI::SUMA_core`/`libSUMA_core` are the documented headless choice.
- CMake and Make build/install/package manifests agree, and both static/shared variants pass their dependency checks.

## Main risks and mitigations

- **Public type source churn:** `SUMA_define.h` embeds many GUI types in broad public structures. Binary layout may change, but all in-tree callers must continue to compile. Mitigate with opaque GUI-private state and compatibility includes in `SUMA_suma.h`.
- **Hidden core-to-GUI calls:** whole-file classification misses calls such as redisplay, popup, event pumping, and viewer attachment. Mitigate with object-level undefined-symbol checks and a forbidden-symbol CI test.
- **Mixed ownership/double free:** object lifecycle is concentrated in files that also render. Move allocation/free pairs together and add sanitizer coverage.
- **Platform linker differences:** a split that works through transitive linkage on macOS may fail with `--no-undefined` on Linux, or vice versa. Test shared and static builds early.
- **False headless status:** merely removing Motif while retaining GL/GLUT/GLX still prevents a truly portable algorithm library. Use the stricter no-graphics completion criterion.
- **Scope growth:** redesigning every SUMA structure at once would be risky. Start with a useful closed core slice, migrate representative consumers, and enlarge it iteratively.

## Resolved design contract

- Compatibility level: source, CLI, input/output, and numerical behavior; not binary ABI.
- Full artifact: `libSUMA`/`AFNI::SUMA`, containing core plus GUI symbols.
- Headless artifact: `libSUMA_core`/`AFNI::SUMA_core`, containing no graphical dependencies.
- Build systems: CMake and legacy Make are both supported.
- Program scope: all `3d*` and `Surf*` programs retain behavior; the directly affected inventory above must run against the headless core.
- Permitted core dependencies: non-graphical AFNI/data/algorithm libraries are allowed.

Implementation should stop and revise this plan only if object/symbol analysis proves that one of the required programs intentionally presents a graphical UI. Otherwise, the program is considered headless and graphical calls in its current dependency path are refactoring work, not justification for leaving it on full `libSUMA`.
