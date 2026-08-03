# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

ECGPACK is a collection of closely related parallel Fortran codes for variational calculations of quantum few-body systems (few-electron atoms and molecules) using all-particle explicitly correlated Gaussian (ECG) basis functions. Developed by the research group of Sergiy Bubin (Physics Department, Nazarbayev University). Parallelism is via MPI. The public repository is <https://github.com/ecgpack/ecgpack> and its default branch is `main`. See `README.md` and directory `doc/` for documentation that includes physics conventions and references.

## Repository layout

Each code lives in its own top-level directory with an identical structure: a `Makefile`, a `src/` subdirectory, and a `.vscode/` subdirectory (nearly identical, user-independent JSON configs for building/debugging the code in VSCode). Most, but not all, codes also include a `sample_input/` subdirectory of worked examples (see below). The codes fall into several groups (some of these groups may overlap):

- **Energy codes (diagonal matrix element codes; basis generation codes)** — `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`, `CG_0S`. These are used for generating ECG basis sets for states of different angular-momentum/parity. They can compute energies, expectation values, particle distributions, and save wave functions.
- **Off-diagonal matrix element codes** — `RG_0S-1P`, `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_1P-2D`, `RG_1P-2P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These codes are used to compute off-diagonal matrix elements between two states of different angular-momentum/parity (and which are thus expanded using different basis sets) for such operators as transition dipole moment, spin–orbit interaction, and non-contact spin–spin interaction.
- **Transition dipole codes** — `RG_0S-1P`, `RG_1P-2D`, `RG_1P-2P`. These codes are used to compute off-diagonal matrix elements of the transition dipole moment operator in the length and velocity gauges.
- **Fine structure coupling codes** — `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These codes are used to compute off-diagonal matrix elements of the spin–orbit interaction and non-contact spin–spin interaction.
- **Real Gaussian codes (Real ECG codes)** — `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`, `RG_0S-1P`, `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_1P-2D`, `RG_1P-2P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These deal with real Gaussian basis sets (real ECGs).
- **Complex Gaussian codes (Complex ECG codes)** — `CG_0S`. These deal with complex ECG basis sets (complex ECGs). Currently there is only one basis ($L=0$), and it is a work in progress that lags behind the real ECG codes. It may not even be compilable.

Most codes ship with a `sample_input/` subdirectory of worked examples — the four real-ECG energy codes (`RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`) plus all nine off-diagonal codes (`RG_0S-1P`, `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_1P-2D`, `RG_1P-2P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`). Only `CG_0S` currently has no sample inputs. See the Running section below.

`RG_0S` is the most complete reference implementation; other codes share much of its source and Makefile. Non-code directories: `doc/` (documentation), `utilities/`, and `bin/` + `jobs/` (created by the user, not in git). The root also holds `README.md` (the repository manual), `AUTHORS.md` (list of contributors), and `.code-workspace` (a VSCode multi-folder workspace grouping all the code directories).

## Building

The canonical entry point is `build.bash` in the root directory — run it with no arguments for full usage. It loops over toolchains/configs/codes/precisions/linalg choices, builds via each code's Makefile, and stores binaries in `bin/<toolchain>/<config>/<CODE>_N<nparticles>_P<precision>_<linalg>` (the `<linalg>` suffix is always present — e.g. `_netlib`, `_mkl`, `_openblas`). The `nparticles` argument is required. The `linalg` argument (see below) takes one of `netlib` (default), `mkl`, `lblas`, `openblas`, `aocl`; for `precision=10`/`16` only `netlib` is built (other values are skipped).

```bash
./build.bash machine=ubuntu-generic toolchain=systemdefault config=release code=RG_0S nparticles=4 precision=8
```

To build a single code directly, invoke its Makefile (this is what `build.bash` calls under the hood):

```bash
cd RG_0S
make release COMPILER=gfortran MACHINE=ubuntu-generic PREC=8 LINALG=openblas EXEFILE=ecg
make debug   COMPILER=gfortran MACHINE=ubuntu-generic PREC=8 LINALG=netlib   EXEFILE=ecg
make clean   # also: cleaner, cleanest, cleanrelease, cleandebug
```

Key build parameters:

- `CONFIG` (`release`/`debug`) — set by the `make release`/`make debug` target. `release` uses `-O3 -march=native`; `debug` enables bounds/uninit/FPE checks. Object and `.mod` files go in `release/` or `debug/`.
- `PREC` — real `kind`: `8` (double/fp64), `10` (extended/fp80, GNU only), `16` (quadruple). Selects which `src/wp_def_<PREC>.f90` is compiled.
- `LINALG` — selects which BLAS/LAPACK implementation to link against (default `netlib`). `netlib` compiles the bundled, lightly modified reference `src/BLAS.f`/`src/LAPACK.f` and adds no extra link flags; `mkl` (Intel MKL — compiler-dependent `-lmkl_*` flags), `lblas` (`-llapack -lblas`), `openblas` (`-lopenblas`), and `aocl` (AMD AOCL `-lflame -lblis …`) instead link an external optimized library and skip the bundled sources. Only `PREC=8` honors the optimized choices; for `PREC=10`/`16` only `LINALG=netlib` is supported (any other value leaves the build unsupported). In the off-diagonal codes `LINALG` is accepted but a no-op (they link no BLAS/LAPACK). The bundled `src/BLAS.f`/`src/LAPACK.f` and the `BARE_OBJS_LPKBLS` object list are compiled only when `LINALG=netlib`.
- `COMPILER` (`gfortran`→`mpif90`, `ifort`→`mpiifort`, `ifx`→`mpiifx`, `nvfortran`→`mpif90`) and `MACHINE` select compiler flags. Supported machines are hardcoded in both `build.bash` and the Makefiles; adding a machine means editing both.

Note: the **number of particles is compiled in**, not a runtime argument. `build.bash` does an in-place `sed` on `src/wp_def_<PREC>.f90` to set `Glob_AllowedNumOfParticles`, builds, then restores the original from `wp_def_temporary.f90`. A binary built for N particles rejects input files with a different particle count.

## Running

Each binary is an MPI program. The energy codes read a single input/output file named **`inout.txt`** from the current working directory (default `Glob_DataFileName` in `globvars.f90`); the file's first (optional) line is `BASIS_TYPE <BASIS>`, the second (required) line is `PARTICLES <n>`. The off-diagonal matrix-element codes instead read two wave-function files (one per state) and do not use `inout.txt`. The transition dipole codes (`RG_0S-1P`, `RG_1P-2D`, `RG_1P-2P`), as well as `RG_0S-2D`, read just the two wave-function files directly, with no extra input file; the remaining fine-structure coupling codes (`RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`) additionally require a small `inp.txt` whose positional argument(s) select the transition/coupling type. The wave-function file names vary by code: most use `wf_state0.txt`/`wf_state1.txt`, but `RG_1P-2D` and `RG_1P-2P` use `wf_state1.txt`/`wf_state2.txt`. Run from a job directory containing the required input file(s):

```bash
mpirun -np <N> /path/to/bin/<toolchain>/<config>/RG_0S_N4_P8_netlib
```

There is no automated test suite; validation is done by running physical test cases and comparing computed energies/expectation values against published reference values. The codes that ship sample inputs (see Repository layout) provide ready-to-run cases under `<CODE>/sample_input/`, documented by `README.md` files giving the description, expected runtime, and the reference energy/values to reproduce. Where those READMEs sit depends on the code: the energy and transition dipole codes put one in each case directory, whereas the fine-structure codes put one in each of the case's three subdirectories (`initial_state_*/`, `final_state_*/`, `matelem/`) and none at the case level. `RG_0S-1P` and `RG_1P-2D` additionally have an overview `README.md` at the top of `sample_input/`. For the energy codes each case is a single directory with an `inout.txt`; recurring kinds are `basis_generation_*` (build an ECG basis from scratch up to a target size), `expected_values_*` (compute expectation values from a saved basis), `densities_*` (compute densities and pair correlation functions, with grid-builder Bash scripts, gnuplot plotting scripts, and sample output plots), and `store_wavefunction_*` (saves both the linear and nonlinear variational parameters of the wave function into a file). The off-diagonal codes' examples instead bundle subdirectories for the two states plus the matrix-element run: the transition dipole cases (`RG_0S-1P`, `RG_1P-2D`, `RG_1P-2P`) are named `transition_dipole_moment_*`, while the fine-structure cases (`RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`) provide `initial_state_*/`, `final_state_*/`, and `matelem/` subdirectories (the `matelem/` directory holds the two wave-function files and, for every fine-structure code except `RG_0S-2D`, the `inp.txt` that selects the coupling type).

## Source architecture

Within each `src/`, the module compile/dependency order (see the Makefile) is:

`wp_def_<PREC>` → `globvars` → `misc`, `linalg`, `spin` → `matelem` → `matform` → `workproc` → `main`

- **`wp_def_<PREC>.f90`** — defines the working-precision kind (`wp`), `Glob_AllowedNumOfParticles`, and the MPI real type. Edited at build time by `build.bash` (see above).
- **`globvars.f90`** — all global state (the `Glob_*` variables: masses, charges, basis, matrices) and physical/numeric constants.
- **`linalg.f90`** — linear-algebra wrappers over BLAS/LAPACK. `BLAS.f`/`LAPACK.f` are bundled, lightly modified netlib reference sources used only when `LINALG=netlib`. `dmng.f` (lightly modified TOMS nonlinear minimizer) and `X1MACH.f90` (machine constants) support the optimizer.
- **`spin.f90`** — spin algebra and permutational-symmetry projection.
- **`matelem.f90`** — matrix elements between individual basis functions; **`matform.f90`** assembles the full Hamiltonian (H) and overlap (S) matrices.
- **`workproc.f90`** — the bulk of the program (often >8000 lines): `ReadIOFile`/`SaveResults` I/O, basis enlargement, optimization cycles, the generalized symmetric eigenvalue solvers (methods `'G'` and `'I'`), expectation values, densities, and swap-file handling.
- **`main.f90`** — initializes MPI, seeds RNGs, then drives a sequence of **BBOP** (Basis Building and Optimization Program) steps read from the input file. Each step is a `case` in the main `select`: `BASIS_ENL`, `OPT_CYCLE`, `FULL_OPT1`, `ELIM_LCFN`, `ELIM_LND1`, `SEPR_LND1`, `SEPR_FLCF`, `EXPC_VALS`, `DENSITIES`, `MOMT_DENS`, `SAVE_FILE`, `SAVE_HSWF`. Adding a calculation mode means adding a case here plus the corresponding routine in `workproc.f90`.

When editing matrix-element or matform code, changes usually must be mirrored across the analogous `RG_*` directories, since the codes are intentionally near-duplicates specialized to different symmetries.

## GPU backend (RG_0S only)

**Read `kb/README.md` before doing any GPU work.** It is the index; it names the current
measurements, the archived-and-wrong ones, and the open questions.

- **`RG_0S/src/gpu_backend.f90`** — CUDA Fortran backend. Exists in **RG_0S only**; `RG_1P`,
  `RG_2P`, `RG_2D` have no GPU path. Two kernels (energy, energy+gradient), one block per
  (k,l) basis pair, 128 threads striding over permutation terms. Both call the *same*
  `MatrixElementsHS_RG_0S` in `matelem.f90`, which carries `attributes(host,device)` and takes
  its physics constants as arguments (device code cannot read `Glob_*` module variables) —
  so **changing that routine's signature breaks the GPU build**.
- Runtime is opt-in and off by default: `ECG_GPU=1` (matrix elements on GPU),
  `ECG_GPU_EIG=1` (cuSOLVER method-G eigensolve; implies `ECG_GPU`), `ECG_GPU_BATCH` (pairs
  per kernel launch, default 16384), `ECG_DETERM=1` (ordered reduction, reproducible),
  `ECG_BENCH=1` (print timings and stop), `ECG_VERIFY=1` (dump H/S/eigenvector),
  `ECG_SEED`, `ECG_GRADCHECK` + `ECG_GC_*`.
- **`build.bash` cannot build a CUDA binary** — it has no `USE_CUDA` knob. GPU builds must call
  the Makefile directly:
  ```bash
  cd RG_0S && make release COMPILER=nvfortran MACHINE=shabyt USE_CUDA=yes PREC=8 \
       LINALG=netlib CUF_MAXREG=128 EXEFILE=ecg
  ```
  `USE_CUDA=yes` requires `COMPILER=nvfortran` and `PREC=8`. `CUDA_ARCH` comes from `MACHINE`
  (shabyt→sm_70, aurora→sm_86).
- The number of particles is compiled in. `make` alone does **not** set it — `sed` it in
  `src/wp_def_8.f90` (`Glob_AllowedNumOfParticles`) or use `build.bash` for CPU builds.

## Benchmarking discipline

Hard-won; ignoring these has produced wrong conclusions in this project before.

1. **Compare against the production build, not a convenient one.** Production runs
   **gfortran + OpenBLAS** (`foss-*` toolchain). Benchmarking the GPU against an nvfortran
   `-Kieee` CPU build overstates it by **~1.75×** — this invalidated six kb documents.
2. **`OMP_NUM_THREADS=1` always.** Multi-threaded OpenBLAS is dramatically *slower* at these
   matrix sizes; parallelise with MPI instead.
3. **Use medians, never means.** The 64-rank CPU arm varies up to 20% within a job and ~28%
   between jobs; the GPU arm is stable to 4 significant figures. A 2-rep mean has already
   produced a 10%-wrong headline.
4. **Time a fixed basis, not a growth run.** `EXPC_VALS` on a stored deck does identical work
   every time; growth trajectories diverge between backends (~1e-9 differences amplify through
   accept/reject steps) and are not comparable run-to-run.
5. **Record the energy in every timing cell.** It is a free correctness gate and has caught
   real problems.
6. **Benchmark the workload you actually run.** Production decks alternate `BASIS_ENL` with
   `OPT_CYCLE`, and the OPT_CYCLE (derivative-ME) half is ~170× the work. Almost all historical
   benchmarking used `BASIS_ENL`-only decks.
7. **A cell with zero work is a FAIL, not an OK.** Gate on `me_calls > 0`, and check that every
   arm of a comparison did the *same* number of calls. Harnesses have twice marked
   `0.000 s over 0 calls` cells `OK` by testing for an empty field instead of a zero one.

## Deck anatomy (the OPT_CYCLE trap)

An `inout.txt` has four `====`-delimited sections: **header**, a single **state line**, the
**BBOP block**, the **K-line history**, then the **basis**. The state line and the last history
line each carry a `CyclesDone` field and **they routinely disagree.**

`main.f90:132-153` guards OPT_CYCLE on the *last history* entry:
`Glob_History(Glob_CurrBasisSize)%CyclesDone < NumCycles`. A harness that reads `CyclesDone`
from the *state* line (section 1) can compute a `NumCycles` that fails the guard — and OPT_CYCLE
is then **silently skipped**: the run completes cleanly, exit 0, zero ME calls. Always take it
from section 3:

```bash
c0=$(awk '/^ *=+ *$/{s++;next} s==3{c=$3} END{print c+0}' "$deck"); nc=$(( c0 + 1 ))
```

This cost three misdiagnoses (blamed the cycle guard, then frozen functions, then "malformed
decks") before anyone read the deck layout and the guard side by side. See
`kb/gradient_path.md` §5.

## Cluster workflow (shabyt)

- Two accounts: `shabyt` (your own) and `shabyt-sb` (Bubin's, holds `~/aidyn/`). **Campaigns
  have been run and left there without ever being copied into `data/`** — always sync results
  back, or they are invisible to everyone including a fresh checkout. See
  `kb/cluster_due_diligence.md`.
- **Never run two build jobs in the same source tree.** Each `sed`s `wp_def_8.f90` and wipes
  `release/`; concurrent jobs corrupt each other and produce spurious build failures.
- **Profiling: NVHPC ships no profiler on shabyt.** Its `profilers/` directory is empty, so
  `which ncu` finds a stub that dies with `ncu-Error-Version 13_0 ... Nsight_Compute is not
  available in this installation`. The working profilers live in the **standalone `CUDA`
  modules**. Load NVHPC for the runtime libraries but call `ncu` by absolute path, matching the
  CUDA version NVHPC was built against:
  ```bash
  module load NVHPC/25.9-CUDA-12.9.1
  NCU=/shared/opt/easybuild/software/CUDA/12.9.1/nsight-compute-2025.2.1/ncu
  ```
  (also present: `CUDA/12.8.0` → `nsight-compute-2025.1.0`, `CUDA/12.6.0` → `2024.3.0`.)
  The same applies to `nsys`.
- **Do not wrap `ncu` in `mpirun`** — it fails with `hwloc_set_cpubind`. Run the binary
  directly (OpenMPI singleton init). Shorten profiled launches with `ECG_GPU_BATCH`, or an
  Oxygen-scale kernel takes far too long to replay.
- GPU performance counters may need elevated permissions (`ERR_NVGPUCTRPERM`). If shabyt
  refuses, use the aurora workstation, where counter access can be enabled locally.
- Do not `module purge` and swap MPI mid-job; it hangs at `MPI_Init`. Set the environment once
  at job start, or use separate jobs.
- `kb/` and `data/` are excluded via `.git/info/exclude` — on disk, not in the shared repo.
