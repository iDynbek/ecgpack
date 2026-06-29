// cuSOLVER backend for the method-G generalized symmetric-definite eigensolve.
// Drop-in GPU replacement for the rank-0 LAPACK DSYGVX calls in workproc.f90
// (EnergyGA / EnergyGAM / EnergyGB).  Solves  H x = lambda S x  (itype 1),
// returning the iwhich-th eigenvalue (1-based, ascending) and, when
// jobz_vec/=0, its eigenvector (normalized so x^T S x = 1, matching DSYGVX).
//
// Uses cusolverDnDsygvdx (divide & conquer, index range), the GPU analogue of
// LAPACK DSYGVX with RANGE=I.  Handle and device buffers are cached across calls
// and grown on demand, so steady-state timing reflects the solver, not allocation
// (mirroring the preallocated Glob_WorkForDSYGVX on the CPU side).
//
// Matrices are Fortran column-major with leading dim ld>=n; cuSOLVER is also
// column-major, so H/S are copied (de-strided to ld=n) with cudaMemcpy2D.

#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <stdio.h>
#include <stdlib.h>

// --- error checking: turn silent OOM / status failures into file:line aborts ---
#define CUDA_CHECK(call) do {                                              \
    cudaError_t _e = (call);                                               \
    if (_e != cudaSuccess) {                                               \
        fprintf(stderr, "CUDA error %s:%d: %s\n",                         \
                __FILE__, __LINE__, cudaGetErrorString(_e));               \
        fflush(stderr); abort();                                          \
    } } while (0)
#define CUSOLVER_CHECK(call) do {                                          \
    cusolverStatus_t _s = (call);                                          \
    if (_s != CUSOLVER_STATUS_SUCCESS) {                                   \
        fprintf(stderr, "cuSOLVER error %s:%d: status %d\n",              \
                __FILE__, __LINE__, (int)_s);                              \
        fflush(stderr); abort();                                          \
    } } while (0)

static cusolverDnHandle_t g_h   = NULL;
static double *g_A = NULL, *g_B = NULL, *g_W = NULL, *g_work = NULL;
static int    *g_info = NULL;
static int     g_cap_n = 0;
static size_t  g_cap_work = 0;

extern "C"
void gpu_dsygvx_(int *jobz_vec, int *n_, double *H, int *ldH_,
                 double *S, int *ldS_, int *iwhich,
                 double *eval, double *Z, int *info_out)
{
    const int n   = *n_;
    const int ldH = *ldH_;
    const int ldS = *ldS_;

    if (g_h == NULL) CUSOLVER_CHECK(cusolverDnCreate(&g_h));
    if (g_info == NULL) CUDA_CHECK(cudaMalloc((void**)&g_info, sizeof(int)));

    if (n > g_cap_n) {
        if (g_A) CUDA_CHECK(cudaFree(g_A));
        if (g_B) CUDA_CHECK(cudaFree(g_B));
        if (g_W) CUDA_CHECK(cudaFree(g_W));
        CUDA_CHECK(cudaMalloc((void**)&g_A, (size_t)n*n*sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&g_B, (size_t)n*n*sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&g_W, (size_t)n*sizeof(double)));
        g_cap_n = n;
    }

    // De-stride copy H,S (host ld -> device n) onto the device, column-major.
    CUDA_CHECK(cudaMemcpy2D(g_A, (size_t)n*sizeof(double), H, (size_t)ldH*sizeof(double),
                            (size_t)n*sizeof(double), (size_t)n, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy2D(g_B, (size_t)n*sizeof(double), S, (size_t)ldS*sizeof(double),
                            (size_t)n*sizeof(double), (size_t)n, cudaMemcpyHostToDevice));

    cusolverEigType_t  itype = CUSOLVER_EIG_TYPE_1;           // H x = lambda S x
    cusolverEigMode_t  jobz  = (*jobz_vec) ? CUSOLVER_EIG_MODE_VECTOR
                                           : CUSOLVER_EIG_MODE_NOVECTOR;
    cusolverEigRange_t range = CUSOLVER_EIG_RANGE_I;
    cublasFillMode_t   uplo  = CUBLAS_FILL_MODE_UPPER;
    const int il = *iwhich, iu = *iwhich;
    int meig = 0, lwork = 0;

    CUSOLVER_CHECK(cusolverDnDsygvdx_bufferSize(g_h, itype, jobz, range, uplo, n,
        g_A, n, g_B, n, 0.0, 0.0, il, iu, &meig, g_W, &lwork));

    if ((size_t)lwork > g_cap_work) {
        if (g_work) CUDA_CHECK(cudaFree(g_work));
        CUDA_CHECK(cudaMalloc((void**)&g_work, (size_t)lwork*sizeof(double)));
        g_cap_work = (size_t)lwork;
    }

    CUSOLVER_CHECK(cusolverDnDsygvdx(g_h, itype, jobz, range, uplo, n,
        g_A, n, g_B, n, 0.0, 0.0, il, iu, &meig, g_W, g_work, lwork, g_info));
    CUDA_CHECK(cudaDeviceSynchronize());

    int info = 0;
    CUDA_CHECK(cudaMemcpy(&info, g_info, sizeof(int), cudaMemcpyDeviceToHost));
    *info_out = info;

    // Selected eigenvalue is the first of the meig returned (ascending).
    CUDA_CHECK(cudaMemcpy(eval, g_W, sizeof(double), cudaMemcpyDeviceToHost));

    // Eigenvector (if requested) is the first column of the overwritten A.
    if (*jobz_vec)
        CUDA_CHECK(cudaMemcpy(Z, g_A, (size_t)n*sizeof(double), cudaMemcpyDeviceToHost));
}

// Free cached device buffers + handle (call from gpu_finalize path at shutdown).
extern "C"
void gpu_eig_finalize_(void)
{
    if (g_A)    { cudaFree(g_A);    g_A = NULL; }
    if (g_B)    { cudaFree(g_B);    g_B = NULL; }
    if (g_W)    { cudaFree(g_W);    g_W = NULL; }
    if (g_work) { cudaFree(g_work); g_work = NULL; }
    if (g_info) { cudaFree(g_info); g_info = NULL; }
    if (g_h)    { cusolverDnDestroy(g_h); g_h = NULL; }
    g_cap_n = 0; g_cap_work = 0;
}
