# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

ECG is a collection of closely related parallel Fortran codes for variational calculations of quantum few-body systems (few-electron atoms and molecules) using all-particle explicitly correlated Gaussian (ECG) basis functions. Developed by the research group of Sergiy Bubin (Physics Department, Nazarbayev University). Parallelism is via MPI. See `README.md` for the physics conventions and references.

## Repository layout

Each code lives in its own top-level directory with an identical structure: a `Makefile` plus a `src/` subdirectory. The codes fall into two groups:

- **Energy/wavefunction codes** — `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P` (real Gaussians for different angular-momentum/parity symmetries) and `CG_0S` (complex Gaussians, L=0; work in progress, lags behind `RG_0S`).
- **Matrix-element codes** — `RG_0S-1P`, `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_1P-2D`, `RG_1P-2P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P` (offdiagonal elements between two bases: transition dipole moments, spin–orbit, spin–spin).

`RG_0S` is the most complete reference implementation; other codes share much of its source and Makefile. Non-code directories: `archive/` (old versions, do not touch), `utilities/`, and `bin/` + `jobs/` (created by the user, not in git).

## Building

The canonical entry point is `build.bash` in the root directory — run it with no arguments for full usage. It loops over toolchains/configs/codes/precisions, builds via each code's Makefile, and stores binaries in `bin/<toolchain>/<config>/<CODE>_N<nparticles>_P<precision>[_optlapack]`. The `nparticles` argument is required.

```bash
./build.bash machine=ubuntu-generic toolchain=systemdefault config=release code=RG_0S nparticles=4 precision=8
```

To build a single code directly, invoke its Makefile (this is what `build.bash` calls under the hood):

```bash
cd RG_0S
make release COMPILER_TYPE=gnu MACHINE=ubuntu-generic PREC=8 OPTLPKBLS=yes EXEFILE=ecg
make debug   COMPILER_TYPE=gnu MACHINE=ubuntu-generic PREC=8 OPTLPKBLS=no  EXEFILE=ecg
make clean   # also: cleaner, cleanest, cleanrelease, cleandebug
```

Key build parameters:

- `CONFIG` (`release`/`debug`) — set by the `make release`/`make debug` target. `release` uses `-O3 -march=native`; `debug` enables bounds/uninit/FPE checks. Object and `.mod` files go in `release/` or `debug/`.
- `PREC` — real `kind`: `8` (double/fp64), `10` (extended/fp80, GNU only), `16` (quadruple). Selects which `src/wp_def_<PREC>.f90` is compiled.
- `OPTLPKBLS` — `yes` links the system optimized LAPACK/BLAS (or MKL for Intel); `no` compiles the bundled reference `src/BLAS.f` and `src/LAPACK.f`. Always forced to `no` for `PREC=10`/`16`.
- `COMPILER_TYPE` (`gnu`→`mpif90`, `intel`→`mpiifort`) and `MACHINE` select compiler flags. Supported machines are hardcoded in both `build.bash` and the Makefiles; adding a machine means editing both.

Note: the **number of particles is compiled in**, not a runtime argument. `build.bash` does an in-place `sed` on `src/wp_def_<PREC>.f90` to set `Glob_MaxAllowedNumOfParticles`, builds, then restores the original from `wp_def_temporary.f90`. A binary built for N particles rejects input files with a different particle count.

## Running

Each binary is an MPI program. It reads a single input/output file named **`inout.txt`** from the current working directory (default `Glob_DataFileName` in `globvars.f90`); the file's first line is `PARTICLES <n>`. Run from a job directory containing `inout.txt`:

```bash
mpirun -np <N> /path/to/bin/.../RG_0S_N4_P8_optlapack
```

There is no automated test suite; validation is done by running physical test cases and comparing computed energies/expectation values against published reference values.

## Source architecture

Within each `src/`, the module compile/dependency order (see the Makefile) is:

`wp_def_<PREC>` → `globvars` → `misc`, `linalg`, `spin` → `matelem` → `matform` → `workproc` → `main`

- **`wp_def_<PREC>.f90`** — defines the working-precision kind (`dprec`), `Glob_MaxAllowedNumOfParticles`, and the MPI real type. Edited at build time by `build.bash` (see above).
- **`globvars.f90`** — all global state (the `Glob_*` variables: masses, charges, basis, matrices) and physical/numeric constants.
- **`linalg.f90`** — linear-algebra wrappers over BLAS/LAPACK. `BLAS.f`/`LAPACK.f` are bundled netlib reference sources used only when `OPTLPKBLS=no`. `dmng.f` (TOMS nonlinear minimizer) and `X1MACH.f90` (machine constants) support the optimizer.
- **`spin.f90`** — spin algebra and permutational-symmetry projection.
- **`matelem.f90`** — matrix elements between individual basis functions; **`matform.f90`** assembles the full Hamiltonian (H) and overlap (S) matrices.
- **`workproc.f90`** — the bulk of the program (often >8000 lines): `ReadIOFile`/`SaveResults` I/O, basis enlargement, optimization cycles, the generalized symmetric eigenvalue solvers (methods `'G'` and `'I'`), expectation values, densities, and swap-file handling.
- **`main.f90`** — initializes MPI, seeds RNGs, then drives a sequence of **BBOP** (Basis Building and Optimization Program) steps read from the input file. Each step is a `case` in the main `select`: `BASIS_ENL`, `OPT_CYCLE`, `FULL_OPT1`, `ELIM_LCFN`, `ELIM_LND1`, `SEPR_LND1`, `SEPR_FLCF`, `EXPC_VALS`, `DENSITIES`, `MOMT_DENS`, `SAVE_FILE`, `SAVE_HSWF`. Adding a calculation mode means adding a case here plus the corresponding routine in `workproc.f90`.

When editing matrix-element or matform code, changes usually must be mirrored across the analogous `RG_*` directories, since the codes are intentionally near-duplicates specialized to different symmetries.
