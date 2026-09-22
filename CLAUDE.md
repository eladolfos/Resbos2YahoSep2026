# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

ResBos2 is a C++11 Monte Carlo physics program that computes fully differential cross-sections
for hadron collider processes (Drell-Yan, W/Z, Higgs, diboson) including QCD soft-gluon
resummation effects. This is research physics code (MSU PhD research), not a typical application —
correctness of the physics/math is paramount, and changes to calculation code should preserve
existing numerical behavior unless a physics change is explicitly intended.

## Build

Out-of-source CMake build. There is already a configured `build/` directory on this machine
(built against HPC-installed LHAPDF/HOPPET under `/home/fuyao3/GCC13/...`), configured with:

```bash
cmake .. -DLHAPDF_ROOT_DIR=/home/fuyao3/GCC13/LHAPDF/ -DHoppet_ROOT_DIR=/home/fuyao3/GCC13/HOPPET/ \
         -DUSE-ROOT=ON -DFITTING=ON -DCMAKE_INSTALL_PREFIX=.
```

(see `build/compile.sh`). To build after configuring:

```bash
cd build
make -jN
make -jN install   # copies binaries/configs to CMAKE_INSTALL_PREFIX
```

Required external dependencies: LHAPDF (PDF sets), HOPPET (Fortran, splitting-function
convolutions — requires `gfortran`), pthreads. `external/` and `CMake/CPM.cmake` pull in
smaller header-only/vendored deps (fmt, yaml-cpp, GTest, CLI11) automatically via CPM.

Key CMake options (`CMake/ProgramOptions.cmake`), all OFF by default except as noted:
- `USE-ROOT` — histogram output via ROOT (currently ON in the local build)
- `USE-YODA` — histogram output via YODA
- `USE-GZIP` — gzip-compressed grid/output I/O
- `USE-MPI` — MPI-parallel grid generation/VEGAS integration
- `USE-OPENMP` — OpenMP parallelism (used by the fitting executables)
- `PYEXT` — build the SWIG Python extension
- `FITTING` — build the non-perturbative fitting executables (`NonPertFit`, `GetFixed`);
  **requires `USE-ROOT=ON`** (enforced with a `FATAL_ERROR` in CMake) and also pulls in
  BAT and Eigen3

Tests (GTest, via CTest) are **only built when `CMAKE_BUILD_TYPE MATCHES DEBUG`**
(see the guard at the bottom of the top-level `CMakeLists.txt`) — a Release configure
(like the existing `build/`) will not build or run tests at all. To build/run tests,
reconfigure a separate debug build directory:

```bash
mkdir build-debug && cd build-debug
cmake .. -DCMAKE_BUILD_TYPE=DEBUG -DLHAPDF_ROOT_DIR=... -DHoppet_ROOT_DIR=...
make -jN
ctest                       # run all test suites
ctest -R testUtilities      # run one suite (testUtilities/testBeam/testProcess/testCalculation/...)
./tests/Utilities/testUtilities --gtest_filter=Grid3D*   # run a single GTest case directly
```

Test binaries are organized to mirror `src/`: `tests/Utilities`, `tests/Beam`, `tests/Process`,
`tests/Calculation`, plus `tests/IntegrationTests` for cross-module tests. Each has its own
`CMakeLists.txt` listing its `*_test.cc` sources.

## Running the code

```bash
build/bin/resbos -i resbos.config -o <output_dir>
```

`resbos.config` (INI-style key/value, parsed by `IO::ConfigParser`/`IO::Settings`) drives
everything about a run: process type, beam/ECM, resummation orders/scheme, non-perturbative
form, kinematic cuts, and histogram mode. See `README.md` for the full list of config keys —
it documents every setting and CLI flag and is the reference for what a given key does.

Fitting the non-perturbative function (requires `-DFITTING=ON`) uses a separate config/data
set: `fitting.config`, `parameters.yml`, `experiments.yml`, `ptData.txt` (copied into the
build dir by CMake when `FITTING=ON`). The `NonPert/` directory holds an older/standalone
fitting driver (`NonPertFit.cxx`, `runNonPert`) kept alongside the CMake-integrated
`src/NonPertFit` version.

## Architecture

### Layered module structure (`src/`, mirrored in `include/ResBos/`)

Static libraries link bottom-up; each is a separate CMake target:

- **`Utilities`** (`ResbosUtilities`) — no internal dependencies. Math/numerics (Vegas
  integrator, adaptive integrators, `Grid1D/2D/3D` interpolation grids, `AdaptiveMap`,
  Gauss-Kronrod/Newton-Cotes/Ogata-quadrature), kinematics primitives (`FourVector`,
  `ThreeVector`, `Lorentz`, `Particle`, `Event`), `Electroweak` parameters,
  `ConfigParser`/`Settings` (input-card parsing), `Histogram`, LHEF event output,
  `ThreadPool`, `Random`, `loguru`-based logging.
- **`Beam`** (`ResbosBeam`, depends on Utilities) — `PDF` (LHAPDF wrapper), `HoppetInterface`
  (splitting-function convolutions via the Fortran HOPPET library), `Convolution`
  (PDF ⊗ splitting-function ⊗ Wilson-coefficient grids), `Beam` (beam/ECM bookkeeping).
- **`Calculation`** (`ResbosCalculation`, depends on Utilities) — the resummation/fixed-order
  math: `Resummation` (CSS/CFG A/B/C/H coefficients), `Asymptotic`, `DeltaSigma`,
  `Perturbative`, `WmA`. Each calculation mode (`XSec` setting) is a subclass registered
  with `CalculationFactory` via the `REGISTER(calculation)` macro
  (see `include/ResBos/Calculation.hh`) — dispatched at runtime by the string in
  `resbos.config`'s `XSec` key.
- **`Process`** (`ResbosProcess`, depends on Utilities) — per-final-state physics:
  `DrellYan`, `Z0`, `A0`, `Wpm`, `Higgs`. Each is a subclass registered with
  `ProcessFactory` via `REGISTER(process)` (`include/ResBos/Process.hh`), dispatched by
  the `Process` config key. A process supplies phase space generation, matrix elements
  (`ME`), non-perturbative form, and `MakeEvent` (fills an `Event` for LHE output).
- **`User`** (`ResbosUser`) — user-customizable cuts/histogram booking
  (`UserCuts.cc`, `UserHistograms.cc`) — the intended place to add analysis-specific cuts
  or histograms without touching core physics code.
- **`NonPertFit`** (only if `FITTING=ON`) — standalone executables `NonPertFit` (fits the
  non-perturbative function to data via BAT) and `GetFixed`, built against all of the above
  libraries plus `ResBos.cc` directly (not linked against the `Resbos` library target).

**To add a new process or calculation type**: subclass `ResBos::Process` /
`ResBos::Calculation`, implement the pure-virtual interface, and add
`REGISTER(YourClassName)` at file scope in the `.cc` — this self-registers into the
factory singleton at static-init time, keyed by the class's `GetName()`. Add the new
`.cc` to the relevant `src/*/CMakeLists.txt` `target_sources`.

### Top-level orchestration

`ResBos::ResBos` (`src/ResBos.cc`, `include/ResBos/ResBos.hh`) is the god object wiring
together `Settings`, `Beam` pair, `Process`, and `Calculation` for one run; `main.cc`
constructs it, then drives grid generation (interpolation grids for PDF convolutions and
the resummed/asymptotic calculation) and finally VEGAS integration over the process phase
space to get the cross-section (`vegas(resbos->GetXSect())`), optionally writing LHE events.

`Enums.hh` centralizes the enums (orders, schemes, convolution types) used across module
boundaries — check there first when threading a new option through multiple layers.

`MCFMInterface` (separate library target in `src/CMakeLists.txt`) exposes ResBos's
asymptotic/resummed pieces for interfacing with the external MCFM NLO generator (see the
`WmA` calculation mode, "used for MCFM interface" per `OrganizationDetails.md`).

### Fortran interop

HOPPET (splitting-function evolution/convolution) and `CT14Pdf.f` are Fortran; `HoppetInterface.cc`
and the CMake `gfortran` link requirement are the bridge points. Be careful with Fortran/C++
calling-convention and array-ordering issues when touching this boundary.

### Design reference

`OrganizationDetails.md` documents the intended class/dependency map (which class depends on
which, and where each is mainly used) — consult it before restructuring dependencies between
modules, since it reflects the intended architecture, not just current-state coupling.
