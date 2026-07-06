# ECG RG_0S — GPU Acceleration Status Report

*Branch `gpu-integration` (tip `97b35d6`). Report compiled 2026-07-06.*

This document summarizes the current state of the GPU-accelerated `RG_0S` code,
the changes introduced on the `gpu-integration` branch, and — in the appendices —
a deep dive into the CUDA kernels and the benchmark results measured to date.

For the operational reference (build flags, runtime knobs, determinism notes) see
[`gpu_acceleration.md`](gpu_acceleration.md); this report is the higher-level
state-of-play plus the performance appendix.

---

## 1. Executive summary

`RG_0S` (real explicitly-correlated Gaussian, L=0 energy code) has an optional
CUDA backend that offloads (a) the **matrix-element build** — the Hamiltonian/
overlap assembly over basis-function pairs — and (b) the **method-G generalized
eigensolve** (via cuSOLVER). Both are opt-in at run time and have no effect on a
CPU build.

Status: **the backend is consolidated, validated, and modular.** It builds and
runs on an NVIDIA V100 (Shabyt) through the full nvfortran + nvcc + cuSOLVER
toolchain, and all GPU-specific code is isolated so the shared Fortran sources
carry only a thin dispatch. Numerically the kernel is exact-to-machine-precision
on a fixed basis; run-to-run it varies at the round-off level (atomicAdd order —
see §A.4).

Headline performance (V100, measured on the current kernels):

- **Matrix elements:** ~**40×** a single CPU core for Carbon, but a flat launch-
  overhead floor means it *loses* for light atoms. Break-even ≈ **Lithium**.
  Against a fully MPI-parallel socket, **1 V100 ≈ one 64-core EPYC node** for the
  ME build.
- **Method-G eigensolve (cuSOLVER):** crossover vs the best CPU at **K ≈ 300**,
  reaching **5.2×** at K = 1000 and still widening — and the CPU eigensolve does
  not parallelize, so this advantage is durable.
- **Production caveat:** the default production path is method **I** (inverse
  iteration), which is *eigensolve-bound* on the serial `GSEPIIS` routine
  (72–80 % of runtime) — a path the current GPU work does **not** yet accelerate.
  That is the highest-value remaining lever (see §5).

---

## 2. Current state of the code

### 2.1 What the GPU backend consists of

| Layer | Files | Compiled |
|---|---|---|
| CUDA matrix-element kernels (energy + gradient) | `src/matelem_cuda.cu` | nvcc, `USE_CUDA` only |
| cuSOLVER method-G eigensolver | `src/eigen_cuda.cu` | nvcc, `USE_CUDA` only |
| `iso_c_binding` interfaces | `src/matelem_gpu.f90` | nvfortran, `USE_CUDA` only |
| **Orchestration** (selection, lifecycle, chunked builds, eigensolve dispatch) | `src/gpu_backend.f90` | nvfortran, `USE_CUDA` only |
| Shared sources (thin dispatch only) | `globvars.f90`, `matform.f90`, `workproc.f90`, `main.f90` | always |

### 2.2 Modular architecture

All GPU-specific Fortran that is **not** a raw C binding lives in
**`gpu_backend.f90`** (compiled only with `-DUSE_CUDA`). It owns run-time
selection, the device-context lifecycle, both chunked matrix-element builds, and
the eigensolver dispatch. It calls back into `matform`'s `StoreHS`/`StoreHSD`
through **procedure-argument callbacks**, so it reuses those storage routines
without a circular module dependency.

Consequently the footprint on the shared sources is near-minimal — only the
call sites that structurally *must* be there:

| Shared file | Lines vs master | What remains |
|---|---|---|
| `globvars.f90` | **0** (byte-identical) | nothing — no GPU state at all |
| `main.f90` | **+6** | one `use` + one collective `call gpu_backend_init()` |
| `matform.f90` | **+18 / −0** | a 5-line early-return dispatch per build routine; **CPU pair-loops byte-identical to master** |
| `workproc.f90` | **~+43** | method-G eigensolve dispatch (`gpu_eig_active()` → cuSOLVER / LAPACK) + a legacy deck-line skip |

Device teardown is automatic: `gpu_init_` registers `atexit(gpu_finalize_)`, so
`main.f90` never calls finalize.

### 2.3 Build and run

```bash
# Build (GPU): CUDA_ARCH is inferred from MACHINE (shabyt->sm_70, aurora->sm_86)
module load NVHPC/25.9-CUDA-12.9.1
make release COMPILER=nvfortran MACHINE=shabyt PREC=8 LINALG=netlib \
     USE_CUDA=yes EXEFILE=ecg_cuda

# Run: GPU is opt-in via env vars; inout.txt stays toolchain-agnostic
ECG_GPU=1                 mpirun -np 1 ./ecg_cuda   # GPU matrix elements
ECG_GPU=1 ECG_GPU_EIG=1   mpirun -np 1 ./ecg_cuda   # + cuSOLVER method-G eigensolve
```

Runtime knobs (all env, all no-ops in a CPU build): `ECG_GPU`, `ECG_GPU_EIG`,
`ECG_GPU_PERSIST` (persistent device buffers, default on), `ECG_GPU_BATCH`
(chunk size, default 16384). GPU builds are `PREC=8` + `LINALG=netlib` only.

---

## 3. Changes introduced on this branch

Five commits, oldest → newest:

1. **Consolidation + docs** (`ab2dc2e`) — merged the CUDA backend onto current
   `master`, reconciled naming/build conventions, ported the cuSOLVER eigensolver
   and CUDA Phases 0/1/2, added `gpu_acceleration.md` + a validation sweep.
2. **Removed the OpenACC prototype** (`f4b94c0`) — the energy-only OpenACC backend
   was a third copy of the physics with no path to parity; archived to branch
   `gpu-openacc-prototype`.
3. **Removed the profiling instrumentation** (`c1cda1a`) — the `system_clock`/
   `ProfAccum` scaffolding was interleaved into hot production routines; archived
   to branch `gpu-profiling`.
4. **UX: env-selected GPU** (`22a0347`) — replaced the `USE_GPU` **deck keyword**
   (which was also written back into every saved deck) with the `ECG_GPU` env var,
   and made `CUDA_ARCH` default from `MACHINE`.
5. **Modularization** (`97b35d6`) — introduced `gpu_backend.f90` and reduced the
   shared-file footprint to the table in §2.2.

### CUDA memory/robustness phases (folded in via commit 1)

- **Phase 0** — error checking (`CUDA_CHECK`/`CUSOLVER_CHECK`) on every CUDA/
  cuSOLVER call.
- **Phase 1** — persistent device context: device buffers are cached and grown
  across calls instead of malloc/free per call (`ECG_GPU_PERSIST`).
- **Phase 2** — batch chunking: the K(K+1)/2 pair triangle is processed in
  fixed-size chunks (`ECG_GPU_BATCH`) so neither host nor device memory scales
  with the full basis — this fixed a large-K out-of-memory wall (the gradient
  buffers were tens of GB at K ≈ 14000).

### Validation status

- **CPU build** compiles/links; He 1Se sample reproduces the reference energy.
- **`-DUSE_CUDA` compile pass** (gfortran) over the whole chain — the abstract-
  interface callbacks compile clean.
- **Shabyt V100 (gn04), job 4361** — full nvfortran + nvcc + cuSOLVER build (arch
  auto-resolved to `sm_70`); ran He 1Se with `ECG_GPU=1` (E = −2.90074) and with
  `ECG_GPU=1 ECG_GPU_EIG=1` (cuSOLVER, E = −2.89656); both to completion with
  clean teardown. The modularization is confirmed end-to-end on real hardware.

---

## 4. What is isolated vs what touches master

The design goal was **"as necessary as possible, not as little as possible"**:
master's existing files carry only code that structurally must live there.

- **Fully isolated** (vanish on a CPU build): `matelem_cuda.cu`, `eigen_cuda.cu`,
  `matelem_gpu.f90`, `gpu_backend.f90`.
- **Necessarily shared** (the dispatch touchpoints): the `matform` early-returns
  to the GPU build, the `workproc` eigensolve branch, the `main` subsystem init.
  Everything else was pushed into `gpu_backend`.

The eigensolve dispatch deliberately stays in `workproc` — it owns the LAPACK
`DSYGVX` workspace, so moving it would couple `gpu_backend` to LAPACK.

---

## 5. Known gaps and next steps

1. **GPU method-I eigensolve ("Phase 3")** — the biggest production lever. The
   production path (method I) is eigensolve-bound on `GSEPIIS` (hand-written
   scalar Fortran: incremental LDLᵀ `LDLTF` + triangular solves `LDLTS` + matvec
   `MTMV`), which is serial (rank-0 only) and does not use BLAS. Porting it to the
   device (or a Fortran-native rewrite with `attributes(device)`) is the durable
   win, and the hardest task.
2. **CUDA Fortran evaluation** — a scoped spike to replace the C++ kernels with
   native CUDA Fortran, which would remove the nvcc/C++/FFI seam and — more
   importantly — lower the barrier to (1). Tracked as a backlog item.
3. **Broaden to other `RG_*` codes** — the ME kernel pattern is mechanical to
   port to the analogous symmetry-specialized codes.
4. **Ship it** — open a PR merging `gpu-integration` onto `master` (conflict-free;
   `origin/master` is a strict ancestor).

---

## Appendix A — CUDA kernel deep dive

The physics: for each ordered basis-function pair (k, l) with k ≥ l, the code
sums the Hamiltonian and overlap matrix elements over the **Y⁺Y permutation
terms** (the symmetry projector), weighted by `YHYCoeff`. The number of terms,
`NumYHYTerms`, grows steeply with the number of identical electrons and is the
dominant cost multiplier for heavy atoms. This is the axis the CUDA backend
parallelizes.

### A.1 Launch geometry — block-per-pair, thread-per-term

```
matelem_batch_energy_kernel<<< npairs, NumYHYTerms >>>   // (energy)
matelem_batch_grad_kernel   <<< npairs, NumYHYTerms, smem >>>   // (energy+gradient)
```

- **`gridDim.x = npairs`** — one CUDA **block per (k, l) pair** in the chunk.
- **`blockDim.x = NumYHYTerms`** (capped at 1024) — one **thread per permutation
  term**.

This mirrors the existing CPU MPI decomposition (which splits over permutation
terms), so the same work partition is expressed on-device. The choice is
deliberate: the per-element linear algebra is tiny (n × n dense with n ≤ 7), far
too small for cuSOLVER-batched routines — a bespoke one-thread-per-element kernel
is the right granularity.

### A.2 Per-element work — `compute_matelem_full` (device function)

Each thread evaluates one (k, l, term) matrix element entirely in registers/local
memory via `__device__ compute_matelem_full(...)`. For an n × n system it builds
`Lk, Ll → Ak, Al`, forms `tAkl`, does a **Cholesky factorization + inverse**,
matrix products and traces for the kinetic term, and a particle-pair loop for the
potential (and, in the gradient kernel, the S/T- and Vkl-gradient slabs). Because
n ≤ 7, the work matrices fit in registers/shared memory; `real*8` throughout
(there is no FP32 path).

### A.3 MPI split preserved inside the kernel

Each thread guards its term with the same modulo test the CPU code uses:

```c
int term = j + 1;
if (((q + term) % NumOfProcs) == ProcID) { ... compute + accumulate ... }
```

So under multi-rank runs each rank's kernel computes only its share of the terms;
the partial `Hout`/`Sout` are then combined with a per-chunk `MPI_ALLREDUCE` on
the host. The device kernel and the host reduction together reproduce the CPU
result.

### A.4 Reduction and determinism — `atomicAdd`

Threads accumulate their weighted contributions into block-shared scalars:

```c
__shared__ double sh_H, sh_S;              // energy kernel
atomicAdd(&sh_H, coeff * Hkl);
atomicAdd(&sh_S, coeff * Skl);
// grad kernel adds dynamic shared slabs sh_Dk[2np], sh_Dl[2np] the same way
```

Thread 0 then writes `sh_H/sh_S` to global `Hout[pair]/Sout[pair]`. The
`atomicAdd` reduction is **not order-deterministic** — floating-point summation
order depends on thread scheduling, so matrix elements vary at the ~1e-15 level
run to run. ECG's ill-conditioned overlap magnifies this to ~1e-8 per eigenvalue,
and basis generation amplifies it to ~1e-6 in a final energy. **This is expected,
not a bug**: correctness is checked on a *fixed* basis (see §B.4), where the
kernel is exact to machine precision, not on bit-identical generation runs.

Note (`cc < 6.0`): FP64 `atomicAdd` in shared memory needs compute capability
≥ 6.0; `gpu_init` warns otherwise. V100 (cc 7.0) and RTX 3060 Ti (cc 8.6) are fine.

### A.5 Gradient kernel specifics

`matelem_batch_grad_kernel` uses **dynamic shared memory** of `2 + 4*np` doubles
(`sh_H, sh_S`, then `Dksum[2np]`, `Dlsum[2np]`), atomic-accumulates the gradient
slabs alongside H/S, and honors a per-pair `grad_l_flag` (whether the l-gradient
is needed). The `(2*np × npairs)` gradient output buffers are the memory wall at
large K — the reason Phase-2 chunking exists.

### A.6 cuSOLVER eigensolver — `eigen_cuda.cu`

`gpu_dsygvx_` wraps **`cusolverDnDsygvdx`** (`itype=1`, i.e. `H x = λ S x`, range
`I`, upper) as a drop-in replacement for the three rank-0 `DSYGVX` calls in
`workproc.f90` (eigenvalue-only 'N', and eigenvector 'V'). The handle and device
buffers are cached and grown across calls, so steady-state timing reflects the
solver rather than allocation (mirroring the preallocated `Glob_WorkForDSYGVX`).
H/S are de-strided onto the device with `cudaMemcpy2D` (column-major, ld → n); the
eigenvector is returned normalized `xᵀSx = 1`, matching DSYGVX. Because the
eigensolve is already rank-0-only, only rank 0 touches the GPU for it.

### A.7 Device memory lifecycle

- Buffers live in a persistent context (Phase 1); `ECG_GPU_PERSIST=0` forces the
  legacy malloc/free-per-call path for A/B testing.
- `gpu_backend` chunks the pair triangle (Phase 2) so a single kernel launch never
  stages more than `ECG_GPU_BATCH` pairs; every rank uses the same cap so the
  per-chunk `MPI_ALLREDUCE` stays collective.
- Teardown (`gpu_matelem_finalize_`, `gpu_eig_finalize_`, `cudaDeviceReset`) runs
  via `atexit`.

---

## Appendix B — Benchmark results

**Provenance & caveat.** The numbers below were measured on Shabyt V100s during
the benchmark study (2026-06-14/15, `~/aidyn/bench`) on the *pre-refactor* code.
The refactor moved orchestration only — **the `.cu` kernels are byte-identical** —
so the performance figures carry over unchanged. All ME timings are the method-G
`ProfAccum` per-call timers; the production method-I path was not instrumented in
that study (its profiling came later — §B.5).

### B.0 The test systems

| Atom | Particles | Pseudoparticles n | `NumYHYTerms` |
|---|---|---|---|
| He | 3 | 2 | 2 |
| Li | 4 | 3 | 6 |
| Be | 5 | 4 | 20 |
| B  | 6 | 5 | 120 |
| C  | 7 | 6 | 720 |

`NumYHYTerms` (permutation terms) grows steeply with identical electrons and is
the per-element cost multiplier — hence GPU wins scale with atomic number.

### B.1 Matrix-element kernel vs one CPU core (aligned window)

Per-call ME build time, identical K = 200→300 window, single rank (2026-06-15,
`aligned`):

| Atom | @ K = 300 CPU (1 core) | @ K = 300 GPU (V100) | Result |
|---|---|---|---|
| He | 0.120 ms | 0.378 ms | **CPU 3.2× faster** |
| Li | 0.437 ms | 0.396 ms | **crossover — GPU 1.1×** |
| C  | 167 ms | 3.95 ms | **GPU 42×** |

The GPU has a flat ~0.4 ms launch-overhead floor, so it only pays off when
`NumYHYTerms × work-per-element` ≫ launch latency. **Break-even ≈ Lithium
(n = 4).** (An earlier "40× for He" figure was an artifact of measuring He at
K = 10; the aligned window corrects it.)

### B.2 Matrix elements: GPU vs a full multicore socket

Carbon ME build @ K = 290 (2026-06-15, `cpuscale`) — the fair comparison, since
the ME build is embarrassingly parallel and the CPU code already MPI-splits over
terms:

| Cores (MPI) | 1 | 2 | 4 | 8 | 16 | 32 | 64 |
|---|---|---|---|---|---|---|---|
| ME time (ms) | 185 | 86.5 | 43.3 | 22.2 | 11.7 | 6.77 | 4.13 |

**1 V100 = 3.86 ms** for the same build. So the "40×" against one core collapses
to **1 V100 ≈ a full 64-core EPYC node** (break-even ≈ 60 cores). Parallel
efficiency is ~100 % to 16 cores, 85 % at 32, 70 % at 64.

### B.3 Method-G eigensolve — cuSOLVER vs CPU (the large-K lever)

Carbon, method-G growth, single rank, GPU matrix elements in all runs (so only the
eigensolver varies). Per-call eigensolve time (2026-06-15, `eigsolve`):

| K | GPU cuSOLVER | OpenBLAS (best CPU, 1 thread) | GPU × vs OpenBLAS | ref LAPACK | GPU × vs ref |
|---|---|---|---|---|---|
| 200 | 3.2 ms | 2.1 ms | 0.7× (CPU wins) | 4.9 ms | 1.6× |
| 300 | 4.9 ms | 5.1 ms | **1.0× crossover** | 15 ms | 3.2× |
| 500 | 8.0 ms | 17.6 ms | 2.2× | 70 ms | 8.6× |
| 700 | 13 ms | 41 ms | 3.1× | 187 ms | 14× |
| 1000 | 20 ms | 106 ms | **5.2×** | 542 ms | **27×** |

Scaling in this range: **GPU ~K^1.19, OpenBLAS ~K^2.48, reference LAPACK
~K^2.94.** The low GPU exponent is overhead/bandwidth-bound below K ≈ 1000 (not
yet compute-saturated O(K³)); the point is the GPU advantage **widens with K**.
Crossover vs the best CPU ≈ K = 300; at K = 1000 the V100 is 5.2× and climbing —
and the CPU eigensolve neither MPI- nor thread-parallelizes, so this win is
durable.

**LAPACK choice matters:** single-thread OpenBLAS is ~3–4× the bundled reference
LAPACK on `DSYGVX` and free (just linking). **Multi-thread OpenBLAS is *slower***
at these K (threading overhead on few-hundred-size matrices) — keep
`OPENBLAS_NUM_THREADS=1` and parallelize via MPI.

### B.4 Correctness — fixed-basis recompute (machine precision)

GPU vs CPU vs published `CURRENT_ENERGY`, single `EXPC_VALS` on a saved basis
(2026-06-15, `xcheck`):

| Atom | K | GPU energy (a.u.) | Agreement |
|---|---|---|---|
| He | 1000 | −2.9033045577287062 | 15–16 figures (machine precision) |
| Li | 500 | −7.4774518226180753 | 15–16 figures |
| Be | 500 (trimmed, n = 5) | −14.665974482418418 | == CPU exactly |

On a clean fixed basis the GPU kernel is exact to machine precision. The ~1e-6
spreads seen during *basis generation* are trajectory divergence from the
atomicAdd round-off (§A.4), not kernel error.

### B.5 Regime and the production caveat

- **Method G (full DSYGVX):** heavy atoms are ME-bound at small K; the GPU ME
  kernel then flips every system into the eigensolve-bound regime past a modest K
  (crossover K\* ≈ He 90 / Li 100 / C 270). This is why the GPU eigensolver
  matters for the G path.
- **Method I (production default):** profiling this session (Carbon, K ≈ 205–210,
  single rank, GPU ME on) showed the inverse-iteration eigensolve **dominates** —
  `BASIS_ENL I` = 80 % eig / 20 % ME, `OPT_CYCLE I` = 72 % eig / 28 % ME. The
  method-I eigensolve is `GSEPIIS` (hand-written scalar Fortran, serial), which the
  current GPU work does **not** accelerate. So on the production path the GPU ME
  kernel speeds up the *minority* of runtime; the durable production win is a GPU
  method-I eigensolve (§5, item 1).

### B.6 One-line takeaways

- GPU matrix elements: worth ~a 64-core node; only for **Li and heavier**.
- GPU method-G eigensolve: **5×** the best CPU at K = 1000, widening, and
  unbeatable by CPU parallelism (which doesn't exist for it).
- The kernel is numerically exact on a fixed basis; run-to-run variation is
  expected round-off, not error.
- The production (method-I) bottleneck is a serial CPU eigensolve the GPU does not
  yet touch — the next real lever.
