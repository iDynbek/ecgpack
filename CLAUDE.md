# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

ECG is a collection of closely related parallel Fortran codes for variational calculations of quantum few-body systems (few-electron atoms and molecules) using all-particle explicitly correlated Gaussian (ECG) basis functions. Developed by the research group of Sergiy Bubin (Physics Department, Nazarbayev University). Parallelism is via MPI. See `README.md` and directory `doc/` for documentation that includes physics conventions and references.

## Repository layout

Each code lives in its own top-level directory with an identical structure: a `Makefile`, a `src/` subdirectory, and a `.vscode/` subdirectory (identical, user-independent JSON configs for building/debugging the code in VSCode). Most, but not all, codes also include a `sample_input/` subdirectory of worked examples (see below). The codes fall into several groups (some of these groups may overlap):

- **Energy codes (diagonal matrix element codes; basis generation codes)** — `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`, `CG_0S`. These are used for generating ECG basis sets for states of different angular-momentum/parity. They can compute energies, expectation values, particle distributions, and save wave functions.
- **Off-diagonal matrix element codes** — `RG_0S-1P`, `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_1P-2D`, `RG_1P-2P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These codes are used to compute off-diagonal matrix elements between two states of different angular-momentum/parity (and which are thus expanded using different basis sets) for such operators as transition dipole moment, spin–orbit interaction, and non-contact spin–spin interaction.
- **Transition dipole codes** — `RG_0S-1P`, `RG_1P-2D`, `RG_1P-2P`. These codes are used to compute off-diagonal matrix elements of the transition dipole moment operator in the length and velocity gauges.
- **Fine structure coupling codes** — `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These codes are used to compute off-diagonal matrix elements of the spin–orbit interaction and non-contact spin–spin interaction.
- **Real Gaussian codes (Real ECG codes)** — `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`, `RG_0S-1P`, `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_1P-2D`, `RG_1P-2P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These deal with real Gaussian basis sets (real ECGs).
- **Complex Gaussian codes (Complex ECG codes)** — `CG_0S`. These deal with complex ECG basis sets (complex ECGs). Currently there is only one basis ($L=0$), and it is a work in progress that lags behind the real ECG codes. It may not even be compilable.

Only some codes ship with a `sample_input/` subdirectory of worked examples — the four energy codes (`RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`) plus `RG_0S-1P` and `RG_1P-2D`. The remaining codes (`CG_0S` and the other off-diagonal codes) currently have no sample inputs. See the Running section below.

`RG_0S` is the most complete reference implementation; other codes share much of its source and Makefile. Non-code directories: `archive/` (old versions, do not touch), `doc/` (documentation), `utilities/`, and `bin/` + `jobs/` (created by the user, not in git). The root also holds `README.md` (the repository manual), `AUTHORS.md` (list of contributors), and `.code-workspace` (a VSCode multi-folder workspace grouping all the code directories).

## Building

The canonical entry point is `build.bash` in the root directory — run it with no arguments for full usage. It loops over toolchains/configs/codes/precisions, builds via each code's Makefile, and stores binaries in `bin/<toolchain>/<config>/<CODE>_N<nparticles>_P<precision>[_optlapack]`. The `nparticles` argument is required.

```bash
./build.bash machine=ubuntu-generic toolchain=systemdefault config=release code=RG_0S nparticles=4 precision=8
```

To build a single code directly, invoke its Makefile (this is what `build.bash` calls under the hood):

```bash
cd RG_0S
make release COMPILER=gfortran MACHINE=ubuntu-generic PREC=8 OPTLPKBLS=yes EXEFILE=ecg
make debug   COMPILER=gfortran MACHINE=ubuntu-generic PREC=8 OPTLPKBLS=no  EXEFILE=ecg
make clean   # also: cleaner, cleanest, cleanrelease, cleandebug
```

Key build parameters:

- `CONFIG` (`release`/`debug`) — set by the `make release`/`make debug` target. `release` uses `-O3 -march=native`; `debug` enables bounds/uninit/FPE checks. Object and `.mod` files go in `release/` or `debug/`.
- `PREC` — real `kind`: `8` (double/fp64), `10` (extended/fp80, GNU only), `16` (quadruple). Selects which `src/wp_def_<PREC>.f90` is compiled.
- `OPTLPKBLS` — `yes` links the system optimized LAPACK/BLAS (or MKL for Intel); `no` compiles the bundled reference `src/BLAS.f` and `src/LAPACK.f`. Always forced to `no` for `PREC=10`/`16`.
- `COMPILER` (`gfortran`→`mpif90`, `ifort`→`mpiifort`) and `MACHINE` select compiler flags. Supported machines are hardcoded in both `build.bash` and the Makefiles; adding a machine means editing both.

Note: the **number of particles is compiled in**, not a runtime argument. `build.bash` does an in-place `sed` on `src/wp_def_<PREC>.f90` to set `Glob_MaxAllowedNumOfParticles`, builds, then restores the original from `wp_def_temporary.f90`. A binary built for N particles rejects input files with a different particle count.

## Running

Each binary is an MPI program. The energy codes read a single input/output file named **`inout.txt`** from the current working directory (default `Glob_DataFileName` in `globvars.f90`); the file's first (optional) line is `BASIS_TYPE <BASIS>`, the second (required) line is `PARTICLES <n>`. The off-diagonal matrix-element codes instead read two wave-function files (one per state, e.g. `wf_state0.txt`/`wf_state1.txt`) and do not need an `inout.txt`. Run from a job directory containing the required input file(s):

```bash
mpirun -np <N> /path/to/bin/<toolchain>/<config>/RG_0S_N4_P8_optlapack
```

There is no automated test suite; validation is done by running physical test cases and comparing computed energies/expectation values against published reference values. The codes that ship sample inputs (see Repository layout) provide ready-to-run cases under `<CODE>/sample_input/`, each in its own directory with an `inout.txt` and a `README.md` describing it, expected runtime, and the reference energy/values to reproduce. Several kinds recur: `basis_generation_*` (build an ECG basis from scratch up to a target size), `expected_values_*` (compute expectation values from a saved basis), `densities_*` (compute densities and pair correlation functions, with grid-builder Bash scripts, gnuplot plotting scripts, and sample output plots), `store_wavefunction_*` (saves both the linear and nonlinear variational parameters of the wave function into a file), and `transition_dipole_moment_*` (compute the transition dipole moment using two wave functions of the initial and final states provided in two separate files).

## Source architecture

Within each `src/`, the module compile/dependency order (see the Makefile) is:

`wp_def_<PREC>` → `globvars` → `misc`, `linalg`, `spin` → `matelem` → `matform` → `workproc` → `main`

- **`wp_def_<PREC>.f90`** — defines the working-precision kind (`wp`), `Glob_MaxAllowedNumOfParticles`, and the MPI real type. Edited at build time by `build.bash` (see above).
- **`globvars.f90`** — all global state (the `Glob_*` variables: masses, charges, basis, matrices) and physical/numeric constants.
- **`linalg.f90`** — linear-algebra wrappers over BLAS/LAPACK. `BLAS.f`/`LAPACK.f` are bundled netlib reference sources used only when `OPTLPKBLS=no`. `dmng.f` (TOMS nonlinear minimizer) and `X1MACH.f90` (machine constants) support the optimizer.
- **`spin.f90`** — spin algebra and permutational-symmetry projection.
- **`matelem.f90`** — matrix elements between individual basis functions; **`matform.f90`** assembles the full Hamiltonian (H) and overlap (S) matrices.
- **`workproc.f90`** — the bulk of the program (often >8000 lines): `ReadIOFile`/`SaveResults` I/O, basis enlargement, optimization cycles, the generalized symmetric eigenvalue solvers (methods `'G'` and `'I'`), expectation values, densities, and swap-file handling.
- **`main.f90`** — initializes MPI, seeds RNGs, then drives a sequence of **BBOP** (Basis Building and Optimization Program) steps read from the input file. Each step is a `case` in the main `select`: `BASIS_ENL`, `OPT_CYCLE`, `FULL_OPT1`, `ELIM_LCFN`, `ELIM_LND1`, `SEPR_LND1`, `SEPR_FLCF`, `EXPC_VALS`, `DENSITIES`, `MOMT_DENS`, `SAVE_FILE`, `SAVE_HSWF`. Adding a calculation mode means adding a case here plus the corresponding routine in `workproc.f90`.

When editing matrix-element or matform code, changes usually must be mirrored across the analogous `RG_*` directories, since the codes are intentionally near-duplicates specialized to different symmetries.
