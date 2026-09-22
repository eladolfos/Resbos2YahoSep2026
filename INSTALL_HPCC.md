# Installing this software on the MSU HPCC (replacing ResBos2)

This documents building and installing this repository on the HPCC, based on your
current `.bashrc`. It replaces the old `ResBos2` install with a separate `ResBos2Yaho`
install, side by side, so the switch is a couple of `.bashrc` edits rather than an
in-place overwrite.

**Config used below** (confirmed choices):
- Install directory: `/mnt/home/lopezels/InstallSources/ResBos2Yaho`
- `FITTING=OFF` — you only need to run resummation cross-sections via the `resbos`
  executable, not fit non-perturbative parameters to data. This means **BAT and Eigen3
  are not needed** (`FITTING` is the only thing in this codebase that requires them).
- `USE-ROOT=ON` — matches your existing `thisroot.sh` sourcing, for histogram output.
- `USE-MPI=ON`, `USE-OPENMP=ON` — for multi-node grid generation/VEGAS integration and
  per-node threading on the HPCC.

## 0. What you already have (from `.bashrc`)

Your `.bashrc` already sets up, and this build will reuse as-is:
- `LHAPDF` at `/mnt/home/lopezels/InstallSources/LHAPDF`
- `HOPPET` at `/mnt/home/lopezels/InstallSources/HOPPET1` (note the `1` suffix)
- `ROOT`, sourced via `thisroot.sh` — **must stay sourced before you configure/build**,
  since CMake's `find_package(ROOT)` reads `$ROOTSYS`/the ROOT CMake config that
  `thisroot.sh` puts on your environment. This is exactly why your existing comment says
  "Debe ir antes de ResBos2" — keep that ordering for the new block too.

`apfel`, the `fiatlux`-related `YAML_CPP_ROOT`/`PKG_CONFIG_PATH` vars, and `MCFM-10.3` are
**not used by this build** — this project fetches and builds its own private copy of
yaml-cpp (via CPM, see below), so your system yaml-cpp install is irrelevant here. Leave
those lines in `.bashrc` untouched; they're for other software.

## 1. Load a toolchain on the HPCC

You need a C++11-capable GCC **with a matching gfortran** (HOPPET is Fortran — Fortran
`.mod` files are not ABI-stable across GCC versions, so the compiler you use here must
match whatever GCC built your `HOPPET1` install), plus CMake ≥ 3.9 and an MPI
implementation (OpenMPI, for `USE-MPI=ON`).

```bash
module purge
module spider GCC          # find the GCC version that matches your HOPPET1 build
module spider OpenMPI       # find an OpenMPI build against that same GCC
module load GCC/<version> OpenMPI/<version>
module load CMake            # or whatever your site calls its CMake module
```

Check `module avail` / `module spider` for exact names on your account — these differ
by HPCC deployment and change over time, so don't hardcode a version without checking.
If you're unsure which GCC built your existing LHAPDF/HOPPET1, check with:

```bash
strings /mnt/home/lopezels/InstallSources/HOPPET1/lib/libhoppet_v1.a | grep -i "GCC:"
```

## 2. Get the source onto the HPCC

The code is now on GitHub at
[eladolfos/Resbos2YahoSep2026](https://github.com/eladolfos/Resbos2YahoSep2026), so clone
it directly on the HPCC rather than copying files over:

```bash
mkdir -p /mnt/home/lopezels/InstallSources
cd /mnt/home/lopezels/InstallSources
git clone https://github.com/eladolfos/Resbos2YahoSep2026.git ResBos2Yaho
cd ResBos2Yaho
```

HPCC login nodes normally have outbound internet, so this should work directly; if a
particular node blocks it, run the clone from a node that does have access (the same one
you'll use for the CPM fetch in the next step).

If the repo is private, either clone over SSH instead (requires an SSH key added to your
GitHub account):

```bash
git clone git@github.com:eladolfos/Resbos2YahoSep2026.git ResBos2Yaho
```

or authenticate over HTTPS with a GitHub personal access token when prompted for a
password.

To update later once you have changes on GitHub:

```bash
cd /mnt/home/lopezels/InstallSources/ResBos2Yaho
git pull
cd build && make -j8 && make -j8 install
```

## 3. Configure

The initial configure needs outbound network access, because CPM (`CMake/CPM.cmake`,
invoked from `external/CMakeLists.txt`) fetches `fmt` and `yaml-cpp` source from GitHub
the first time you run `cmake`. Run this step on a login node (or any node with internet),
**not** inside a compute job that has no network.

```bash
cd /mnt/home/lopezels/InstallSources/ResBos2Yaho
mkdir build && cd build

cmake .. \
    -DLHAPDF_ROOT_DIR=/mnt/home/lopezels/InstallSources/LHAPDF \
    -DHoppet_ROOT_DIR=/mnt/home/lopezels/InstallSources/HOPPET1 \
    -DUSE-ROOT=ON \
    -DUSE-MPI=ON \
    -DUSE-OPENMP=ON \
    -DFITTING=OFF \
    -DCMAKE_INSTALL_PREFIX=/mnt/home/lopezels/InstallSources/ResBos2Yaho
```

Watch the CMake output for `Looking for LHAPDF... found` / `Looking for Hoppet... found`
/ `Building with ROOT histograms` / `Building ResBos with MPI` / `...with OpenMP` — if any
say "not found", double check the loaded modules and the `_ROOT_DIR` paths above.

## 4. Build

Compiling is CPU-heavy (this is a physics Monte Carlo code with heavy templated math),
so do this in an interactive job or as a batch job rather than on a login node:

```bash
salloc --time=01:00:00 --cpus-per-task=8 --mem=16G
cd /mnt/home/lopezels/InstallSources/ResBos2Yaho/build
make -j8
```

(No network needed at this point — CPM already fetched everything during configure.)

## 5. Install

```bash
make -j8 install
```

This populates `/mnt/home/lopezels/InstallSources/ResBos2Yaho/{bin,lib}` — the same
`bin`/`lib` layout your old `ResBos2` block in `.bashrc` already expects, just under the
new directory name.

## 6. Update `.bashrc`

Add a new block for `ResBos2Yaho`, keeping it **after** the `thisroot.sh` line, same as
the old block:

```bash
# Rutas de ResBos2Yaho (replaces ResBos2)
export PATH=/mnt/home/lopezels/InstallSources/ResBos2Yaho/bin:$PATH
export LD_LIBRARY_PATH=/mnt/home/lopezels/InstallSources/ResBos2Yaho/lib:$LD_LIBRARY_PATH
```

You can leave the old `ResBos2` block in place while you verify the new build works, then
remove it (or comment it out) once you've switched over — the two `bin`/`lib` directories
are independent, so having both on `PATH`/`LD_LIBRARY_PATH` briefly is harmless as long as
executable/library names don't collide (they don't: both install a `resbos` binary, so
whichever comes first on `PATH` wins — put `ResBos2Yaho` first, or remove the old block,
once you're ready to switch).

If you enabled `USE-MPI=ON`/`USE-OPENMP=ON`, also load the same `OpenMPI` module in
`.bashrc` (or in your job scripts) that you built against, since the `resbos` binary is
dynamically linked against it. Your existing `OMP_NUM_THREADS=28`/`OMP_STACKSIZE=2G`
settings apply as-is to OpenMP-parallel parts of this build.

## 7. PDF sets

`resbos.config` (the shipped default) uses `PDF = CT14nnlo` via LHAPDF. Make sure that set
is installed for your LHAPDF instance before running:

```bash
lhapdf install CT14nnlo
```

(swap in whatever PDF sets your actual runs need — `lhapdf list` shows what's already
installed).

## 8. Verify

Serial smoke test:

```bash
cd <some scratch/output dir>
resbos -i /mnt/home/lopezels/InstallSources/ResBos2Yaho/resbos.config -o .
```

MPI smoke test (since `USE-MPI=ON`):

```bash
srun -n 4 resbos -i /mnt/home/lopezels/InstallSources/ResBos2Yaho/resbos.config -o .
```

A successful run prints the ResBos splash, logs progress to `resbos.log`
(or `resbos_<rank>.log` under MPI), and ends with a `Total Cross-Section: ... +/- ...`
summary block.

## If you later need the fitting executables (BAT)

If you end up needing `NonPertFit`/`GetFixed` after all, that requires installing BAT
(and Eigen3) on the HPCC, then reconfiguring with `-DFITTING=ON -DUSE-ROOT=ON` and setting
`BAT_ROOT_DIR` (see `CMake/FindBAT.cmake` for the exact search variable) alongside the
LHAPDF/Hoppet paths above — everything else in this document stays the same.
