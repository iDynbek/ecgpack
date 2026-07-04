# GPU Acceleration in RG_0S — Changes, Optimizations, and Validation

> Status: consolidated on branch `gpu-integration`. Applies to **`RG_0S`** only
> (the other `RG_*` codes have no GPU path yet). GPU work is opt-in at runtime;
> a GPU-enabled binary still runs a normal CPU calculation when the GPU is not
> requested.

This document describes every GPU-related change and optimization introduced on
the `gpu-integration` branch, how to build and run them, the runtime controls,
the validation experiments and their results, and the known gaps.

---

## 1. What the branch consolidates

`gpu-integration` = current `master` merged with the CUDA/OpenACC GPU backend,
plus the GPU eigensolver and three CUDA memory/robustness phases, all reconciled
to master's conventions:

| Component | Source files | Runtime switch |
|---|---|---|
| CUDA matrix-element backend (energy + gradient) | `src/matelem_cuda.cu`, `src/matelem_gpu.f90`, `matform.f90` | `USE_GPU` deck line |
| OpenACC matrix-element backend (energy, prototype) | `src/matelem_acc.f90` | `USE_GPU` (OpenACC build) |
| cuSOLVER generalized eigensolver (method **G**) | `src/eigen_cuda.cu`, `matelem_gpu.f90`, `workproc.f90` | `ECG_GPU_EIG=1` |
| Phase 0 — CUDA error checking | `matelem_cuda.cu`, `eigen_cuda.cu` | always on |
| Phase 1 — persistent device context | `matelem_cuda.cu` | `ECG_GPU_PERSIST` (default on) |
| Phase 2 — matrix-element batch chunking | `matform.f90` | `ECG_GPU_BATCH` (default 16384) |
| Method-I profiling instrumentation | `workproc.f90`, `main.f90` | always on (near-zero cost) |

The matrix-element build is MPI-parallel (per-rank partial sums + `MPI_ALLREDUCE`);
the eigensolve runs on rank 0.

---

## 2. Building

The GPU host toolchain is **NVIDIA HPC SDK (nvfortran) + nvcc + HPC-X (nvfortran-built OpenMPI)**.
The number of particles is compiled in (`Glob_AllowedNumOfParticles` in `src/wp_def_<PREC>.f90`).

### CPU build (no GPU)
```bash
make release COMPILER=gfortran MACHINE=ubuntu-generic PREC=8 LINALG=openblas EXEFILE=ecg
```
`matform.f90` now always carries `#ifdef USE_CUDA`, so every compiler is built
with preprocessing: gfortran `-cpp`, ifort/ifx `-fpp`, nvfortran `-Mpreprocess`.

### CUDA GPU build
```bash
# load NVHPC (>=25.3; see the 24.9 caveat below) — its bundled HPC-X mpif90 wraps nvfortran
module load NVHPC/25.9-CUDA-12.9.1
make release COMPILER=nvfortran MACHINE=shabyt PREC=8 LINALG=netlib \
     USE_CUDA=yes CUDA_ARCH=sm_70 EXEFILE=ecg_cuda      # sm_70 V100 / sm_80 A100 / sm_86 RTX 3060 Ti
```
`USE_CUDA=yes` compiles `matelem_cuda.cu` + `eigen_cuda.cu` with nvcc, adds
`-DUSE_CUDA`, and links `-cuda -c++libs -cudalib=cusolver,cublas`. GPU builds are
`PREC=8` only and use `LINALG=netlib` (the bundled reference BLAS/LAPACK compiled
by nvfortran, to avoid linking a system BLAS against nvfortran).

### OpenACC GPU build (mutually exclusive with CUDA)
```bash
make release COMPILER=nvfortran MACHINE=shabyt PREC=8 LINALG=netlib \
     USE_OPENACC=yes OPENACC_ARCH=cc70 EXEFILE=ecg_acc
```

### Toolchain caveats
- **Use NVHPC ≥ 25.3.** nvfortran **24.9 ICEs** compiling `linalg.f90`
  (`Lowering Error: bad ast optype`). 25.9 and 26.3 compile it fine.
- **MPI must be the nvfortran-built one.** nvhpc bundles it as **HPC-X**
  (`comm_libs/hpcx`, an OpenMPI build); `module load NVHPC/25.9…` puts its `mpif90`
  in `PATH`. A GCC-built OpenMPI produces an `mpi.mod` nvfortran cannot read.
- The nvfortran `-O3` build is **slow** (~15 min; it grinds on `linalg.f90`).

---

## 3. Runtime controls

| Control | Where | Default | Effect |
|---|---|---|---|
| `USE_GPU` | deck line, immediately after `GENERATOR_PARAM` | off | Route the matrix-element build to the GPU (CUDA or OpenACC). Read positionally by the parser. |
| `ECG_GPU_EIG` | env var | `0` | `1` = solve method-**G** eigenproblems on the GPU (cuSOLVER) instead of LAPACK `DSYGVX`. |
| `ECG_GPU_PERSIST` | env var | `1` | `1` = persistent/grown device buffers (Phase 1); `0` = legacy malloc/free per call (for A/B). |
| `ECG_GPU_PROFILE` | env var | `0` | `1` = print accumulated GPU time / call counts at finalize. |
| `ECG_GPU_BATCH` | env var | `16384` | Max pairs per GPU matrix-element chunk (Phase 2). Read on rank 0, broadcast to all ranks. |

For a 2-GPU node, run `mpirun -np 2` with `--gres=gpu:2`; device = `rank % nGPUs`.
Use `mpirun --bind-to none` on shabyt to avoid an HPC-X cpu-binding error.

The program also prints a **PROFILING SUMMARY** at the end: total time and call
counts for the matrix-element build vs the eigensolve, and their percentage split.

---

## 4. Changes and optimizations in detail

### 4.1 CUDA matrix-element backend
`matelem_cuda.cu` contains the device kernels for the energy and energy+gradient
matrix elements; `matelem_gpu.f90` is the `iso_c_binding` interface; `matform.f90`
dispatches to them when `Glob_UseGPU` is set. Each basis-function pair is a CUDA
block; threads accumulate the contributions of the Y⁺Y operator terms. This is the
dominant speedup for heavier atoms (≈40× for Carbon in earlier benchmarks).

### 4.2 OpenACC backend
`matelem_acc.f90` is a single-source OpenACC port of the energy path (prototype),
built with `USE_OPENACC=yes`. Mutually exclusive with the CUDA backend.

### 4.3 cuSOLVER GPU eigensolver (method G)
`eigen_cuda.cu` wraps `cusolverDnDsygvdx` (divide-and-conquer, index range) — the
GPU analogue of LAPACK `DSYGVX`. Enabled with `ECG_GPU_EIG=1`; the three method-G
`DSYGVX` call sites in `workproc.f90` (`EnergyGA`/`EnergyGB`) dispatch to it. Only
method **G** is covered (method I uses inverse iteration — see §6). Buffers are
cached and grown on demand.

### 4.4 Phase 0 — error checking
`CUDA_CHECK` / `CUSOLVER_CHECK` macros wrap every allocation, copy, and solve, plus
`cudaGetLastError` + `cudaDeviceSynchronize` after kernels. Turns silent OOM /
launch failures into `file:line` aborts with the requested byte count. No behaviour
change.

### 4.5 Phase 1 — persistent device context
Replaces the per-call `cudaMalloc`/upload/launch/`cudaFree` pattern with cached,
grow-on-demand device buffers. Run-invariant inputs (Y⁺Y operator, masses, charges)
are uploaded once. Eliminates the ~10–12 serializing alloc/free round-trips per call
that dominated latency for the small batches the optimizer issues. Legacy path kept
under `ECG_GPU_PERSIST=0` for A/B comparison.

### 4.6 Phase 2 — batch chunking (large-K OOM fix)
Previously both GPU matrix-element paths staged the **entire** K(K+1)/2 pair triangle
in one allocation, so a full rebuild at large K blew out host and device memory (the
gradient buffers, sized `npt2 × npairs`, were ~tens of GB at K≈14000). Both
`ComputeMatElem` (energy) and `ComputeMatElemAndDeriv` (gradient) now loop over the
pairs in fixed-size chunks (`ECG_GPU_BATCH`, default 16384), allocating buffers sized
to the chunk. A single flush site (chunk full **or** last pair) keeps the per-chunk
`MPI_ALLREDUCE` collective; the cap is broadcast from rank 0 so every rank chunks
identically. Rank 0 logs the pair/chunk counts when chunking engages. With Phase 1's
grow-on-demand device buffers, **host and device memory are now bounded by the chunk
size at any K**. Chunking is numerically exact: each pair's element is independent and
the per-element `ALLREDUCE` sum is unchanged by grouping.

### 4.7 Method-I profiling instrumentation
The profiling timers previously wrapped only the method-G routines. `EnergyIA`/`IAM`/`IB`
now time their `ComputeMatElem[AndDeriv]` (as ME) and `GSEPIIS` (as eigensolve) calls,
so a method-I run reports a real ME-vs-eigensolve split.

---

## 5. Validation results

Reproduce with `utilities/gpu_validation_sweep.sh` (one GPU job, all configs).

### 5.1 Correctness

Consolidated sweep (Carbon, **fixed** 201-function basis via `EXPC_VALS G 201`, so
every config solves the identical eigenproblem; CPU matrix elements are exact and are
the ground truth). Single V100, `NVHPC/25.9`.

| Config | Validates | Energy (hartree) | Δ vs CPU ref |
|---|---|---|---|
| `cpu_ref` | CPU ME + CPU `DSYGVX` (exact) | −37.50158509184552 | — |
| `gpu_me_a` | GPU ME (noise run 1) | −37.50158504004103 | −5.2e-8 |
| `gpu_me_b` | GPU ME (noise run 2) | −37.50158503239668 | −5.9e-8 |
| `gpu_eig` | GPU cuSOLVER eigensolver | −37.50158507578148 | −1.6e-8 |
| `persist_off` | Phase 1 off (legacy malloc/free) | −37.50158500123170 | −9.1e-8 |
| `batch_64` | Phase 2 chunking (318 chunks) | −37.50158508660349 | −0.5e-8 |
| `rank2` | 2 ranks / 2 GPUs | −37.50158504317909 | −4.9e-8 |

**All configs agree to 8 significant figures** (all within ~9e-8, i.e. the `atomicAdd`
eigenvalue-noise magnitude — see §6; the two GPU-ME runs agree to 7.6e-9). This
confirms: the GPU matrix-element build is correct (vs exact CPU), the cuSOLVER
eigensolver is correct, Phase 1 (persistent context) is numerically neutral, Phase 2
(chunking) is exact, and 2-rank/2-GPU is correct.

Additional single-purpose checks run during development:
- GPU eigensolver vs LAPACK `DSYGVX` on an identical basis: agree to **1.4e-8**.
- Phase 2 chunked (`ECG_GPU_BATCH=64`, 318 chunks) vs single batch: within the
  single-vs-single noise floor (2.5e-8 vs 2.8e-8) → chunking introduces no error.

### 5.2 Performance / profiling

Method-I ME-vs-eigensolve split (Carbon K≈205–210, single-rank V100, GPU ME on):

| Workload | Matrix-element (GPU) | Eigensolve (`GSEPIIS`, CPU) | Eigensolve share |
|---|---|---|---|
| `BASIS_ENL I` (generation) | 24.4 s / 5478 calls | 100.1 s / 5478 calls | **80.4 %** |
| `OPT_CYCLE I` (optimization) | 66.5 s / 9841 calls | 169.0 s / 9871 calls | **71.8 %** |

The method-I production path is **eigensolve-bound** even with the GPU ME kernel
active, and the eigensolve is serial (rank-0 only) — see §7.

---

## 6. Numerical determinism (important when comparing runs)

The CUDA kernel sums the Y⁺Y operator terms with `atomicAdd`, whose floating-point
order depends on thread scheduling. Consequences:

- Matrix elements vary at the ~1e-15 (round-off) level run to run.
- ECG's ill-conditioned overlap matrix magnifies this to **~1e-8 per eigenvalue**
  (two identical fixed-basis GPU runs differ this much).
- Basis **generation** amplifies it further to **~1e-6** in the final energy.

So GPU-vs-CPU, 1-rank-vs-2-rank, chunked-vs-unchunked, and persist-on-vs-off
**generation** energies are expected to differ by ~1e-6 — this is *not* a bug. To
check correctness, compare a **fixed-basis** single eigenvalue (agrees to ~1e-8) or
compare against published reference energies to ~8 figures. Do not expect
bit-identical generation energies.

---

## 7. Known gaps and future work

- **Method-I GPU eigensolve (Phase 3).** Profiling (§5.2) shows the method-I
  eigensolve (`GSEPIIS`, inverse iteration) is 72–80 % of production runtime and does
  not parallelize (serial, rank-0 only) — the Amdahl scaling wall. A GPU port would be
  the highest-value next lever, but it is **not** a cuSOLVER wrapper (cuSOLVER is
  method-G's full-decomposition approach): it means porting the incremental LDLᵀ
  factorization (`LDLTF`), triangular solves (`LDLTS`), and matvec (`MTMV`) to device
  kernels and keeping H/S/factors resident. `GSEPIIS` is factorization-bound scalar
  Fortran (not BLAS), so a cheaper first step is to vectorize/block those CPU routines.
- **Phase 3 (resident H/S + in-place eigensolve)** — prerequisite for the large-K
  single-GPU pipeline; collides with the multi-rank partial-sum model (single-rank first).
- **Other codes.** GPU support exists only in `RG_0S`.
- **OpenACC backend** is an energy-only prototype.
