// CUDA backend for MatrixElements, RG_0S.
// Provides batched energy-only and energy+gradient kernels matching
// ComputeMatElem / ComputeMatElemAndDeriv in matform.f90.
//
// All matrices are row-major, 0-indexed inside device code.
// Fortran arrays (NonlinParam, YHYMatr, MassMatrix, P) are column-major;
// indexing accounts for this where noted.

#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// --- error checking: turn silent OOM / launch failures into file:line aborts ---
#define CUDA_CHECK(call) do {                                              \
    cudaError_t _e = (call);                                               \
    if (_e != cudaSuccess) {                                               \
        fprintf(stderr, "CUDA error %s:%d: %s\n",                         \
                __FILE__, __LINE__, cudaGetErrorString(_e));               \
        fflush(stderr); abort();                                          \
    } } while (0)

#define NN 7                    // max pseudoparticles
#define TWO_OVER_SQRTPI 1.1283791670955126

// ---------------------------------------------------------------------------
// ScaledChargeProd inline (mirrors Fortran function in matelem.f90)
// ---------------------------------------------------------------------------
#define SCP(q1,q2,as,rs,rp,rm) \
    ((q1)*(q2) < 0.0 ? (q1)*(q2)*(as) : \
     ((q1)>0.0 && (q2)>0.0 ? (q1)*(q2)*(rs)*(rp) : (q1)*(q2)*(rs)*(rm)))

// ---------------------------------------------------------------------------
// Device function: full MatrixElements for one permutation term.
// Computes energy (always) and optionally gradients (grad_k / grad_l).
//
// Inputs/conventions:
//   Lk_v, Ll_v : vech of Lk, Ll  (lower-triangular, col-major packing)
//   P_col       : n×n permutation matrix, Fortran column-major flat array
//                 P(row,col) = P_col[(col-1)*n + (row-1)] in 0-indexed = P_col[c*n+r]
//   massmat     : n×n M matrix, Fortran column-major: M(r,c) = massmat[c*n+r]
//   charge[n]   : pseudoparticle charges q_1..q_n
//   charge0     : reference particle charge q_0
//
// Outputs:
//   *Hkl_out, *Skl_out : energy & overlap for this j-term
//   Dk_out[2*np]       : [dH/dvechLk, dS/dvechLk]  (if grad_k)
//   Dl_out[2*np]       : [dH/dvechLl, dS/dvechLl]  (if grad_l)
// ---------------------------------------------------------------------------
__device__ void compute_matelem_full(
    const double * __restrict__ Lk_v,
    const double * __restrict__ Ll_v,
    const double * __restrict__ P_col,
    int n, int np,
    double piraised3n2,
    const double * __restrict__ massmat,
    const double * __restrict__ charge,
    double charge0,
    double attr_scale, double rep_scale, double rep_plus, double rep_minus,
    int grad_k, int grad_l,
    double *Hkl_out, double *Skl_out,
    double * __restrict__ Dk_out,
    double * __restrict__ Dl_out
)
{
    // -----------------------------------------------------------------------
    // Local matrices (row-major, stack-allocated)
    // -----------------------------------------------------------------------
    double Lk[NN][NN], Ll[NN][NN];
    double Ak[NN][NN], tAl[NN][NN], tAkl[NN][NN];
    double W1[NN][NN], W2[NN][NN];
    double inv_tAkl[NN][NN], inv_ttAkl[NN][NN];
    double inv_tAkltAlM[NN][NN];
    double F[NN][NN], G[NN][NN];
    double tQ[NN][NN], ttQ[NN][NN];
    double tr_inv_tAklJij32[NN][NN];
    double dHkldvechLk[NN*(NN+1)/2], dSkldvechLk[NN*(NN+1)/2];
    double dHkldvechLl[NN*(NN+1)/2], dSkldvechLl[NN*(NN+1)/2];
    double WVcLk[NN*(NN+1)/2], WVcLl[NN*(NN+1)/2];
    double det_tAkl, Skl, Tkl, Vkl;
    double temp1, temp2, temp3, temp4, temp5;
    int i, j, k, kk, kkk, indx;

    // ===== ENERGY PATH (mirrors matelem.f90 lines 68–242) =====

    // Build Lk, Ll from vech (lower-triangular, col-major packing)
    indx = 0;
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            Lk[i][j] = 0.0;  Lk[j][i] = Lk_v[indx];
            Ll[i][j] = 0.0;  Ll[j][i] = Ll_v[indx];
            indx++;
        }
    }

    // Ak = Lk*Lk^T,  tAl = Ll*Ll^T  (before permutation)
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            double sAk = 0.0, sAl = 0.0;
            for (k = 0; k <= i; k++) {
                sAk += Lk[i][k] * Lk[j][k];
                sAl += Ll[i][k] * Ll[j][k];
            }
            Ak[i][j] = Ak[j][i] = sAk;
            tAl[i][j] = tAl[j][i] = sAl;
        }
    }

    // tAl = P^T * Al * P,  tAkl = Ak + tAl
    // P_col[c*n+r] = P(r+1,c+1) in Fortran (column-major)
    // Step 1: W1 = P^T * tAl → W1[j][i] = sum_k P_col[j*n+k] * tAl[k][i]
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            temp1 = 0.0;
            for (k = 0; k < n; k++) temp1 += P_col[j*n+k] * tAl[k][i];
            W1[j][i] = temp1;
        }
    // Step 2: tAl = W1 * P → tAl[i][j] = sum_k W1[i][k] * P_col[j*n+k]
    for (i = 0; i < n; i++)
        for (j = i; j < n; j++) {
            temp1 = 0.0;
            for (k = 0; k < n; k++) temp1 += W1[i][k] * P_col[j*n+k];
            tAl[i][j] = tAl[j][i] = temp1;
            tAkl[i][j] = tAkl[j][i] = Ak[i][j] + temp1;
        }

    // Cholesky of tAkl → lower triangle of W1
    det_tAkl = 1.0;
    for (i = 0; i < n; i++)
        for (j = i; j < n; j++) {
            temp1 = tAkl[i][j];
            for (k = 0; k < i; k++) temp1 -= W1[i][k] * W1[j][k];
            if (i == j) { W1[i][i] = sqrt(temp1); det_tAkl *= temp1; }
            else         { W1[j][i] = temp1 / W1[i][i]; W1[i][j] = 0.0; }
        }

    // Invert Cholesky factor in-place
    for (i = 0; i < n; i++) {
        W1[i][i] = 1.0 / W1[i][i];
        for (j = i+1; j < n; j++) {
            temp1 = 0.0;
            for (k = i; k < j; k++) temp1 -= W1[j][k] * W1[k][i];
            W1[j][i] = temp1 / W1[j][j];
        }
    }

    // inv_tAkl = W1^T * W1
    for (i = 0; i < n; i++)
        for (j = i; j < n; j++) {
            temp1 = 0.0;
            for (k = j; k < n; k++) temp1 += W1[k][i] * W1[k][j];
            inv_tAkl[i][j] = inv_tAkl[j][i] = temp1;
        }

    // Overlap
    Skl = piraised3n2 / (det_tAkl * sqrt(det_tAkl));

    // W2 = inv_tAkl * tAl
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            temp1 = 0.0;
            for (k = 0; k < n; k++) temp1 += inv_tAkl[j][k] * tAl[k][i];
            W2[j][i] = temp1;
        }

    // inv_tAkltAlM = W2 * M  (M col-major: M(r,c) = massmat[c*n+r])
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            temp1 = 0.0;
            for (k = 0; k < n; k++) temp1 += W2[j][k] * massmat[i*n+k];
            inv_tAkltAlM[j][i] = temp1;
        }

    // Kinetic: Tkl = 6 * Skl * tr(inv_tAkltAlM * Ak)
    Tkl = 0.0;
    for (i = 0; i < n; i++) {
        temp1 = 0.0;
        for (k = 0; k < n; k++) temp1 += inv_tAkltAlM[i][k] * Ak[k][i];
        Tkl += temp1;
    }
    Tkl *= 6.0 * Skl;

    // Potential: Vkl = (2/sqrt(pi))*Skl * sum_pairs scp/sqrt(tr[C^-1 J_ij])
    temp1 = TWO_OVER_SQRTPI * Skl;
    Vkl = 0.0;
    for (i = 0; i < n; i++) {
        temp3 = inv_tAkl[i][i]; temp4 = sqrt(temp3);
        tr_inv_tAklJij32[i][i] = 1.0 / (temp4 * temp3);
        Vkl += SCP(charge[i], charge0, attr_scale, rep_scale, rep_plus, rep_minus) * (temp1/temp4);
    }
    for (i = 0; i < n; i++)
        for (j = i+1; j < n; j++) {
            temp3 = inv_tAkl[i][i] + inv_tAkl[j][j] - 2.0*inv_tAkl[j][i];
            temp4 = sqrt(temp3);
            tr_inv_tAklJij32[j][i] = 1.0 / (temp4 * temp3);
            Vkl += SCP(charge[i], charge[j], attr_scale, rep_scale, rep_plus, rep_minus) * (temp1/temp4);
        }

    *Hkl_out = Tkl + Vkl;
    *Skl_out = Skl;

    if (!grad_k && !grad_l) return;

    // ===== S/T-GRADIENT PATH (matelem.f90 lines 244–483) =====

    double temp_3Skl = -3.0 * Skl;

    if (grad_k) {
        // W1 = inv_tAkl * Lk  (Lk lower-triangular: Lk[k][j]=0 for k<j)
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                temp1 = 0.0;
                for (k = j; k < n; k++) temp1 += inv_tAkl[i][k] * Lk[k][j];
                W1[i][j] = temp1;
            }
        // dSkldvechLk[indx] = -3*Skl * W1[j][i]  for i<=j
        indx = 0;
        for (i = 0; i < n; i++)
            for (j = i; j < n; j++)
                dSkldvechLk[indx++] = W1[j][i] * temp_3Skl;

        // F = inv_tAkltAlM * tAl * inv_tAkl
        // Step 1: W1 = inv_tAkltAlM * tAl
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                temp1 = 0.0;
                for (k = 0; k < n; k++) temp1 += inv_tAkltAlM[i][k] * tAl[k][j];
                W1[i][j] = temp1;
            }
        // Step 2: F = W1 * inv_tAkl (symmetric)
        for (i = 0; i < n; i++)
            for (j = i; j < n; j++) {
                temp1 = 0.0;
                for (k = 0; k < n; k++) temp1 += W1[i][k] * inv_tAkl[k][j];
                F[i][j] = F[j][i] = temp1;
            }
        // W1 = F * Lk (lower triangle of Lk)
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                temp1 = 0.0;
                for (k = j; k < n; k++) temp1 += F[i][k] * Lk[k][j];
                W1[i][j] = temp1;
            }
        // dHkldvechLk
        double t1 = Tkl/Skl, t2 = 12.0*Skl;
        indx = 0;
        for (i = 0; i < n; i++)
            for (j = i; j < n; j++, indx++)
                dHkldvechLk[indx] = t1*dSkldvechLk[indx] + t2*W1[j][i];
    }

    if (grad_l) {
        // inv_ttAkl = P * inv_tAkl * P^T
        // Step 1: W1 = P * inv_tAkl → W1[i][j] = sum_k P_col[k*n+i] * inv_tAkl[k][j]
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                temp1 = 0.0;
                for (k = 0; k < n; k++) temp1 += P_col[k*n+i] * inv_tAkl[k][j];
                W1[i][j] = temp1;
            }
        // Step 2: inv_ttAkl = W1 * P^T → inv_ttAkl[i][j] = sum_k W1[i][k] * P_col[k*n+j]
        for (i = 0; i < n; i++)
            for (j = i; j < n; j++) {
                temp1 = 0.0;
                for (k = 0; k < n; k++) temp1 += W1[i][k] * P_col[k*n+j];
                inv_ttAkl[i][j] = inv_ttAkl[j][i] = temp1;
            }
        // W1 = inv_ttAkl * Ll (Ll lower-triangular)
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                temp1 = 0.0;
                for (k = j; k < n; k++) temp1 += inv_ttAkl[i][k] * Ll[k][j];
                W1[i][j] = temp1;
            }
        // dSkldvechLl
        indx = 0;
        for (i = 0; i < n; i++)
            for (j = i; j < n; j++)
                dSkldvechLl[indx++] = W1[j][i] * temp_3Skl;

        // G = P * inv_tAkl * Ak * M * Ak * inv_tAkl * P^T
        // Step 1: W1 = Ak * M  (M col-major)
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                temp1 = 0.0;
                for (k = 0; k < n; k++) temp1 += Ak[i][k] * massmat[j*n+k];
                W1[i][j] = temp1;
            }
        // Step 2: W2 = W1 * Ak = Ak*M*Ak (symmetric)
        for (i = 0; i < n; i++)
            for (j = i; j < n; j++) {
                temp1 = 0.0;
                for (k = 0; k < n; k++) temp1 += W1[i][k] * Ak[k][j];
                W2[i][j] = W2[j][i] = temp1;
            }
        // Step 3: W1 = inv_tAkl * W2
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                temp1 = 0.0;
                for (k = 0; k < n; k++) temp1 += inv_tAkl[i][k] * W2[k][j];
                W1[i][j] = temp1;
            }
        // Step 4: W2 = W1 * inv_tAkl (symmetric)
        for (i = 0; i < n; i++)
            for (j = i; j < n; j++) {
                temp1 = 0.0;
                for (k = 0; k < n; k++) temp1 += W1[i][k] * inv_tAkl[k][j];
                W2[i][j] = W2[j][i] = temp1;
            }
        // Step 5: W1 = P * W2 → W1[i][j] = sum_k P_col[k*n+i] * W2[k][j]
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                temp1 = 0.0;
                for (k = 0; k < n; k++) temp1 += P_col[k*n+i] * W2[k][j];
                W1[i][j] = temp1;
            }
        // Step 6: G = W1 * P^T (symmetric)
        for (i = 0; i < n; i++)
            for (j = i; j < n; j++) {
                temp1 = 0.0;
                for (k = 0; k < n; k++) temp1 += W1[i][k] * P_col[k*n+j];
                G[i][j] = G[j][i] = temp1;
            }
        // W1 = G * Ll (lower triangle of Ll)
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                temp1 = 0.0;
                for (k = j; k < n; k++) temp1 += G[i][k] * Ll[k][j];
                W1[i][j] = temp1;
            }
        // dHkldvechLl
        double t1 = Tkl/Skl, t2 = 12.0*Skl;
        indx = 0;
        for (i = 0; i < n; i++)
            for (j = i; j < n; j++, indx++)
                dHkldvechLl[indx] = t1*dSkldvechLl[indx] + t2*W1[j][i];
    }

    // ===== Vkl-GRADIENT PATH (matelem.f90 lines 489–603) =====

    if (grad_k) for (k = 0; k < np; k++) WVcLk[k] = 0.0;
    if (grad_l) for (k = 0; k < np; k++) WVcLl[k] = 0.0;

    // Double loop over particle pairs (i,j) with j<=i (Fortran: j from 1 to i)
    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++) {
            // Build tQ = inv_tAkl * J_ij * inv_tAkl
            if (i == j) {
                // J_ii = e_i * e_i^T  →  tQ[k][kk] = inv_tAkl[k][i]*inv_tAkl[i][kk]
                for (k = 0; k < n; k++) {
                    tQ[k][k] = inv_tAkl[k][i] * inv_tAkl[i][k];
                    for (kk = k+1; kk < n; kk++) {
                        double v = inv_tAkl[k][i] * inv_tAkl[i][kk];
                        tQ[k][kk] = tQ[kk][k] = v;
                    }
                }
            } else {
                // J_ij = (e_i - e_j)*(e_i - e_j)^T
                for (k = 0; k < n; k++) {
                    tQ[k][k] = inv_tAkl[k][i]*inv_tAkl[i][k]
                              + inv_tAkl[k][j]*inv_tAkl[j][k]
                              - inv_tAkl[k][i]*inv_tAkl[j][k]
                              - inv_tAkl[k][j]*inv_tAkl[i][k];
                    for (kk = k+1; kk < n; kk++) {
                        double v = inv_tAkl[k][i]*inv_tAkl[i][kk]
                                 + inv_tAkl[k][j]*inv_tAkl[j][kk]
                                 - inv_tAkl[k][i]*inv_tAkl[j][kk]
                                 - inv_tAkl[k][j]*inv_tAkl[i][kk];
                        tQ[k][kk] = tQ[kk][k] = v;
                    }
                }
            }

            if (grad_k) {
                // W1 = tQ * Lk (lower triangle of Lk → k_inner >= kk)
                for (k = 0; k < n; k++)
                    for (kk = 0; kk <= k; kk++) {
                        temp2 = 0.0;
                        for (kkk = kk; kkk < n; kkk++) temp2 += tQ[k][kkk] * Lk[kkk][kk];
                        W1[k][kk] = temp2;
                    }
            }
            if (grad_l) {
                // ttQ = P * tQ * P^T
                // Step 1: W2 = P * tQ → W2[k][kk] = sum_kkk P_col[kkk*n+k] * tQ[kkk][kk]
                for (k = 0; k < n; k++)
                    for (kk = 0; kk < n; kk++) {
                        temp2 = 0.0;
                        for (kkk = 0; kkk < n; kkk++) temp2 += P_col[kkk*n+k] * tQ[kkk][kk];
                        W2[k][kk] = temp2;
                    }
                // Step 2: ttQ = W2 * P^T (symmetric)
                for (k = 0; k < n; k++)
                    for (kk = k; kk < n; kk++) {
                        temp2 = 0.0;
                        for (kkk = 0; kkk < n; kkk++) temp2 += W2[k][kkk] * P_col[kkk*n+kk];
                        ttQ[k][kk] = ttQ[kk][k] = temp2;
                    }
                // W2 = ttQ * Ll (lower triangle of Ll → kkk >= kk)
                for (k = 0; k < n; k++)
                    for (kk = 0; kk <= k; kk++) {
                        temp2 = 0.0;
                        for (kkk = kk; kkk < n; kkk++) temp2 += ttQ[k][kkk] * Ll[kkk][kk];
                        W2[k][kk] = temp2;
                    }
            }

            // Charge scaling factor for this pair
            double cpair = (i == j) ? SCP(charge0, charge[i], attr_scale, rep_scale, rep_plus, rep_minus)
                                    : SCP(charge[i], charge[j], attr_scale, rep_scale, rep_plus, rep_minus);
            temp1 = cpair * tr_inv_tAklJij32[i < j ? j : i][i < j ? i : j];
            // (tr_inv_tAklJij32 stored with larger index first in our convention)

            if (grad_k) {
                // Accumulate W1 lower-triangle into WVcLk (vech order: i=0..n-1, j=0..i)
                indx = 0;
                for (k = 0; k < n; k++) {
                    WVcLk[indx++] += temp1 * W1[k][k];
                    for (kk = k+1; kk < n; kk++)
                        WVcLk[indx++] += temp1 * W1[kk][k];
                }
            }
            if (grad_l) {
                indx = 0;
                for (k = 0; k < n; k++) {
                    WVcLl[indx++] += temp1 * W2[k][k];
                    for (kk = k+1; kk < n; kk++)
                        WVcLl[indx++] += temp1 * W2[kk][k];
                }
            }
        }
    }

    // Final gradient: combine Vkl part with H/S gradients
    double vscale = TWO_OVER_SQRTPI * Skl;
    double vklSkl = Vkl / Skl;
    if (grad_k) {
        for (k = 0; k < np; k++)
            dHkldvechLk[k] += vklSkl * dSkldvechLk[k] + vscale * WVcLk[k];
    }
    if (grad_l) {
        for (k = 0; k < np; k++)
            dHkldvechLl[k] += vklSkl * dSkldvechLl[k] + vscale * WVcLl[k];
    }

    // Pack outputs: Dk = [dH/dvechLk, dS/dvechLk],  Dl = [dH/dvechLl, dS/dvechLl]
    if (grad_k) {
        for (k = 0; k < np; k++) { Dk_out[k] = dHkldvechLk[k]; Dk_out[np+k] = dSkldvechLk[k]; }
    }
    if (grad_l) {
        for (k = 0; k < np; k++) { Dl_out[k] = dHkldvechLl[k]; Dl_out[np+k] = dSkldvechLl[k]; }
    }
}

// ---------------------------------------------------------------------------
// Kernel: batched energy-only.
// gridDim.x = npairs; blockDim.x = NumYHYTerms (capped at 1024).
// Block pair_idx handles one (k,l) pair; threads cover permutation terms.
// ---------------------------------------------------------------------------
__global__ void matelem_batch_energy_kernel(
    const double * __restrict__ NonlinParam, // (np, Kmax) col-major: func k at NonlinParam+k*np
    const int    * __restrict__ k_list,      // 1-indexed k for each pair
    const int    * __restrict__ l_list,      // 1-indexed l for each pair
    const double * __restrict__ YHYMatr,     // n*n*T, Fortran col-major 3D slices
    const double * __restrict__ YHYCoeff,
    int NumYHYTerms, int n, int np,
    int NumOfProcs, int ProcID,
    double piraised3n2,
    const double * __restrict__ massmat,
    const double * __restrict__ charge,
    double charge0,
    double attr_scale, double rep_scale, double rep_plus, double rep_minus,
    double * __restrict__ Hout,
    double * __restrict__ Sout
)
{
    __shared__ double sh_H, sh_S;
    int pair_idx = blockIdx.x;

    if (threadIdx.x == 0) { sh_H = 0.0; sh_S = 0.0; }
    __syncthreads();

    int k0 = k_list[pair_idx] - 1;  // 0-indexed basis function index
    int l0 = l_list[pair_idx] - 1;
    int flat_i = pair_idx + 1;      // 1-indexed flat pair counter (matches Fortran i)
    int q = (flat_i - 1) * NumYHYTerms - 1;

    const double *Lk_v = NonlinParam + k0 * np;
    const double *Ll_v = NonlinParam + l0 * np;

    int j = threadIdx.x;
    if (j < NumYHYTerms) {
        int term = j + 1;
        if (((q + term) % NumOfProcs) == ProcID) {
            const double *P_j = YHYMatr + j * n * n;
            double Hkl, Skl;
            double dummy_Dk[NN*(NN+1)], dummy_Dl[NN*(NN+1)];  // unused for energy-only
            compute_matelem_full(Lk_v, Ll_v, P_j, n, np, piraised3n2,
                                  massmat, charge, charge0,
                                  attr_scale, rep_scale, rep_plus, rep_minus,
                                  0, 0,   // grad_k=0, grad_l=0
                                  &Hkl, &Skl, dummy_Dk, dummy_Dl);
            double coeff = YHYCoeff[j];
            atomicAdd(&sh_H, coeff * Hkl);
            atomicAdd(&sh_S, coeff * Skl);
        }
    }
    __syncthreads();

    if (threadIdx.x == 0) {
        Hout[pair_idx] = sh_H;
        Sout[pair_idx] = sh_S;
    }
}

// ---------------------------------------------------------------------------
// Kernel: batched energy + gradients.
// Dk/Dl arrays are (2*np, npairs) Fortran column-major (i.e., Dk[pair_idx*2*np + ...]).
// ---------------------------------------------------------------------------
__global__ void matelem_batch_grad_kernel(
    const double * __restrict__ NonlinParam,
    const int    * __restrict__ k_list,
    const int    * __restrict__ l_list,
    const int    * __restrict__ grad_l_flag, // 1 if grad_l for this pair, 0 otherwise
    const double * __restrict__ YHYMatr,
    const double * __restrict__ YHYCoeff,
    int NumYHYTerms, int n, int np,
    int NumOfProcs, int ProcID,
    double piraised3n2,
    const double * __restrict__ massmat,
    const double * __restrict__ charge,
    double charge0,
    double attr_scale, double rep_scale, double rep_plus, double rep_minus,
    double * __restrict__ Hout,
    double * __restrict__ Sout,
    double * __restrict__ Dkout,   // (2*np, npairs) Fortran col-major
    double * __restrict__ Dlout    // (2*np, npairs) Fortran col-major
)
{
    extern __shared__ double sh_buf[];  // [2 + 4*np] doubles: H,S,Dksum[2np],Dlsum[2np]
    double *sh_H    = sh_buf;
    double *sh_S    = sh_buf + 1;
    double *sh_Dk   = sh_buf + 2;
    double *sh_Dl   = sh_buf + 2 + 2*np;

    int pair_idx = blockIdx.x;
    int do_grad_l = grad_l_flag[pair_idx];

    // Zero shared memory
    if (threadIdx.x == 0) {
        *sh_H = 0.0; *sh_S = 0.0;
    }
    for (int q = threadIdx.x; q < 2*np; q += blockDim.x) {
        sh_Dk[q] = 0.0;
        sh_Dl[q] = 0.0;
    }
    __syncthreads();

    int k0 = k_list[pair_idx] - 1;
    int l0 = l_list[pair_idx] - 1;
    int flat_i = pair_idx + 1;
    int q_base = (flat_i - 1) * NumYHYTerms - 1;

    const double *Lk_v = NonlinParam + k0 * np;
    const double *Ll_v = NonlinParam + l0 * np;

    int j = threadIdx.x;
    if (j < NumYHYTerms) {
        int term = j + 1;
        if (((q_base + term) % NumOfProcs) == ProcID) {
            const double *P_j = YHYMatr + j * n * n;
            double Hkl, Skl;
            double Dk[NN*(NN+1)], Dl[NN*(NN+1)];
            compute_matelem_full(Lk_v, Ll_v, P_j, n, np, piraised3n2,
                                  massmat, charge, charge0,
                                  attr_scale, rep_scale, rep_plus, rep_minus,
                                  1, do_grad_l,
                                  &Hkl, &Skl, Dk, Dl);
            double coeff = YHYCoeff[j];
            atomicAdd(sh_H, coeff * Hkl);
            atomicAdd(sh_S, coeff * Skl);
            for (int q = 0; q < 2*np; q++) {
                atomicAdd(&sh_Dk[q], coeff * Dk[q]);
                if (do_grad_l) atomicAdd(&sh_Dl[q], coeff * Dl[q]);
            }
        }
    }
    __syncthreads();

    if (threadIdx.x == 0) {
        Hout[pair_idx] = *sh_H;
        Sout[pair_idx] = *sh_S;
        double *Dkp = Dkout + pair_idx * 2*np;
        double *Dlp = Dlout + pair_idx * 2*np;
        for (int q = 0; q < 2*np; q++) {
            Dkp[q] = sh_Dk[q];
            if (do_grad_l) Dlp[q] = sh_Dl[q];
        }
    }
}

// ---------------------------------------------------------------------------
// Persistent device context.
//
// The matrix-element entry points are called once per gradient/energy
// evaluation inside the nonlinear optimizer (workproc.f90 -> dmng), i.e.
// hundreds–thousands of times per optimization window. The original code
// cudaMalloc'd ~12 buffers, re-uploaded *all* inputs (including data that is
// constant for the whole run), launched, then cudaFree'd everything on every
// call. cudaMalloc/cudaFree are synchronous, device-serializing ops, so that
// per-call churn dominated latency for the small batches the optimizer issues.
//
// Here device buffers are cached and grown on demand (mirroring eigen_cuda.cu /
// the preallocated Glob_WorkForDSYGVX on the CPU side), and the run-constant
// inputs (YHYMatr, YHYCoeff, massmat, charge) are uploaded only once.
//
// The legacy per-call path is preserved and selectable at runtime with
//   ECG_GPU_PERSIST=0   (default: 1 = persistent)
// so the two implementations can be A/B compared from the same binary. Set
//   ECG_GPU_PROFILE=1   to print accumulated time/call counts at finalize.
// ---------------------------------------------------------------------------
typedef struct {
    double *d_NLP, *d_YHYMatr, *d_YHYCoeff, *d_massmat, *d_charge;
    double *d_H, *d_S, *d_Dk, *d_Dl;
    int    *d_klist, *d_llist, *d_grfl;
    size_t  c_NLP, c_YHYMatr, c_YHYCoeff, c_massmat, c_charge;
    size_t  c_H, c_S, c_Dk, c_Dl;
    size_t  c_klist, c_llist, c_grfl;
    int     constants_uploaded;     // YHYMatr/YHYCoeff/massmat/charge resident?
} MatElemCtx;

static MatElemCtx g_me = {0};

// Accumulated wall-clock profiling (host timer brackets full sync'd GPU work).
static double  g_prof_energy = 0.0, g_prof_grad = 0.0;
static long    g_calls_energy = 0,  g_calls_grad = 0;

// ECG_GPU_PERSIST: 1 (default) = cache+grow buffers; 0 = legacy malloc/free per call.
static int me_persist(void) {
    static int m = -1;
    if (m < 0) { const char *e = getenv("ECG_GPU_PERSIST"); m = (e && e[0]=='0') ? 0 : 1; }
    return m;
}
static int me_profile(void) {
    static int m = -1;
    if (m < 0) { const char *e = getenv("ECG_GPU_PROFILE"); m = (e && e[0]=='1') ? 1 : 0; }
    return m;
}
static double me_now(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1.0e-9*(double)ts.tv_nsec;
}

// Grow a cached buffer to >= `need` elements; reallocs only when capacity is
// exceeded. In legacy mode capacities are reset to 0 after each call (see
// me_release_buffers), so this allocates fresh every time, reproducing the
// original behavior exactly.
static void me_grow_d(double **p, size_t *cap, size_t need) {
    if (*cap >= need) return;
    if (*p) CUDA_CHECK(cudaFree(*p));
    CUDA_CHECK(cudaMalloc((void**)p, need*sizeof(double)));
    *cap = need;
}
static void me_grow_i(int **p, size_t *cap, size_t need) {
    if (*cap >= need) return;
    if (*p) CUDA_CHECK(cudaFree(*p));
    CUDA_CHECK(cudaMalloc((void**)p, need*sizeof(int)));
    *cap = need;
}

// Free all cached device buffers and reset the context.
static void me_release_buffers(void) {
    double *dp[] = { g_me.d_NLP, g_me.d_YHYMatr, g_me.d_YHYCoeff, g_me.d_massmat,
                     g_me.d_charge, g_me.d_H, g_me.d_S, g_me.d_Dk, g_me.d_Dl };
    for (unsigned i = 0; i < sizeof(dp)/sizeof(dp[0]); i++) if (dp[i]) cudaFree(dp[i]);
    int *ip[] = { g_me.d_klist, g_me.d_llist, g_me.d_grfl };
    for (unsigned i = 0; i < sizeof(ip)/sizeof(ip[0]); i++) if (ip[i]) cudaFree(ip[i]);
    MatElemCtx z = {0};
    g_me = z;
}

extern "C" {

void gpu_matelem_finalize_(void) {
    if (me_profile()) {
        fprintf(stderr,
            "gpu_matelem: energy %ld calls, %.4f s (%.1f us/call); "
            "grad %ld calls, %.4f s (%.1f us/call) [persist=%d]\n",
            g_calls_energy, g_prof_energy,
            g_calls_energy ? 1.0e6*g_prof_energy/g_calls_energy : 0.0,
            g_calls_grad, g_prof_grad,
            g_calls_grad ? 1.0e6*g_prof_grad/g_calls_grad : 0.0,
            me_persist());
        fflush(stderr);
    }
    me_release_buffers();
}

void gpu_init_(int *local_rank)
{
    /* Register device teardown to run automatically at process exit, so the
       Fortran side never has to call gpu_finalize_ explicitly (keeps main.f90
       free of GPU lifecycle code). Guarded against repeated gpu_init_ calls. */
    static int finalize_registered = 0;
    void gpu_finalize_(void);
    if (!finalize_registered) { atexit(gpu_finalize_); finalize_registered = 1; }
    int ndev = 0;
    cudaGetDeviceCount(&ndev);
    if (ndev == 0) { fprintf(stderr, "gpu_init: no CUDA devices\n"); return; }
    int dev = (*local_rank) % ndev;
    cudaSetDevice(dev);
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, dev);
    printf("gpu_init: rank %d → device %d (%s, %.1f GB, cc %d.%d)\n",
           *local_rank, dev, prop.name,
           prop.totalGlobalMem/1.0e9, prop.major, prop.minor);
    if (prop.major < 6)
        printf("gpu_init: WARNING: cc < 6.0 — FP64 atomicAdd in shared memory not supported\n");
    else if (prop.major == 6 && prop.minor == 0)
        printf("gpu_init: NOTE: consumer GPU — FP64 rate ~1/32 of FP32\n");
}

void gpu_finalize_(void) {
    void gpu_eig_finalize_(void);
    gpu_matelem_finalize_();
    gpu_eig_finalize_();
    cudaDeviceReset();
}

// ---- Batched energy-only (replaces ComputeMatElem inner j-loop) ----
// NonlinParam: (np, Nmax) Fortran col-major
// k_list[npairs], l_list[npairs]: 1-indexed
// Hout[npairs], Sout[npairs]: partial sums (this rank's MPI contribution)
void gpu_compute_matelem_batch_(
    const double *NonlinParam, const int *Nmax_in,
    const int *k_list, const int *l_list, const int *npairs_in,
    const double *YHYMatr, const double *YHYCoeff, const int *NumYHYTerms_in,
    const int *n_in, const int *np_in,
    const int *NumOfProcs_in, const int *ProcID_in,
    const double *piraised3n2,
    const double *massmat, const double *charge, const double *charge0,
    const double *attr_scale, const double *rep_scale,
    const double *rep_plus,   const double *rep_minus,
    double *Hout, double *Sout
)
{
    int npairs = *npairs_in, Nmax = *Nmax_in;
    int nn = *n_in, nnp = *np_in, T = *NumYHYTerms_in;
    int persist = me_persist(), profile = me_profile();
    double t0 = profile ? me_now() : 0.0;

    // Ensure capacity (allocates only on growth in persistent mode).
    me_grow_d(&g_me.d_NLP,     &g_me.c_NLP,     (size_t)Nmax*nnp);
    me_grow_d(&g_me.d_YHYMatr, &g_me.c_YHYMatr, (size_t)nn*nn*T);
    me_grow_d(&g_me.d_YHYCoeff,&g_me.c_YHYCoeff,(size_t)T);
    me_grow_d(&g_me.d_massmat, &g_me.c_massmat, (size_t)nn*nn);
    me_grow_d(&g_me.d_charge,  &g_me.c_charge,  (size_t)nn);
    me_grow_d(&g_me.d_H,       &g_me.c_H,       (size_t)npairs);
    me_grow_d(&g_me.d_S,       &g_me.c_S,       (size_t)npairs);
    me_grow_i(&g_me.d_klist,   &g_me.c_klist,   (size_t)npairs);
    me_grow_i(&g_me.d_llist,   &g_me.c_llist,   (size_t)npairs);

    // Per-call inputs (change every call): always uploaded.
    CUDA_CHECK(cudaMemcpy(g_me.d_NLP,   NonlinParam, (size_t)Nmax*nnp*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(g_me.d_klist, k_list,      (size_t)npairs*sizeof(int),      cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(g_me.d_llist, l_list,      (size_t)npairs*sizeof(int),      cudaMemcpyHostToDevice));

    // Run-constant inputs: uploaded once (every call in legacy mode).
    if (!g_me.constants_uploaded) {
        CUDA_CHECK(cudaMemcpy(g_me.d_YHYMatr,  YHYMatr,  (size_t)nn*nn*T*sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(g_me.d_YHYCoeff, YHYCoeff, (size_t)T*sizeof(double),       cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(g_me.d_massmat,  massmat,  (size_t)nn*nn*sizeof(double),   cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(g_me.d_charge,   charge,   (size_t)nn*sizeof(double),      cudaMemcpyHostToDevice));
        g_me.constants_uploaded = 1;
    }

    int threads = min(T, 1024);
    matelem_batch_energy_kernel<<<npairs, threads>>>(
        g_me.d_NLP, g_me.d_klist, g_me.d_llist, g_me.d_YHYMatr, g_me.d_YHYCoeff,
        T, nn, nnp, *NumOfProcs_in, *ProcID_in,
        *piraised3n2, g_me.d_massmat, g_me.d_charge, *charge0,
        *attr_scale, *rep_scale, *rep_plus, *rep_minus,
        g_me.d_H, g_me.d_S
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(Hout, g_me.d_H, npairs*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Sout, g_me.d_S, npairs*sizeof(double), cudaMemcpyDeviceToHost));

    if (!persist) me_release_buffers();
    if (profile) { g_prof_energy += me_now() - t0; g_calls_energy++; }
}

// ---- Batched energy+gradient (replaces ComputeMatElemAndDeriv inner j-loop) ----
// grad_l_list[npairs]: 1 where l > nfru AND l != k, else 0
// Dkout[npairs * 2*np], Dlout[npairs * 2*np]: partial gradient sums
void gpu_compute_matelem_and_deriv_batch_(
    const double *NonlinParam, const int *Nmax_in,
    const int *k_list, const int *l_list,
    const int *grad_l_list,
    const int *npairs_in,
    const double *YHYMatr, const double *YHYCoeff, const int *NumYHYTerms_in,
    const int *n_in, const int *np_in,
    const int *NumOfProcs_in, const int *ProcID_in,
    const double *piraised3n2,
    const double *massmat, const double *charge, const double *charge0,
    const double *attr_scale, const double *rep_scale,
    const double *rep_plus,   const double *rep_minus,
    double *Hout, double *Sout,
    double *Dkout, double *Dlout
)
{
    int npairs = *npairs_in, Nmax = *Nmax_in;
    int nn = *n_in, nnp = *np_in, T = *NumYHYTerms_in;
    int npt2 = 2 * nnp;
    int persist = me_persist(), profile = me_profile();
    double t0 = profile ? me_now() : 0.0;

    // Ensure capacity (allocates only on growth in persistent mode).
    me_grow_d(&g_me.d_NLP,     &g_me.c_NLP,     (size_t)Nmax*nnp);
    me_grow_d(&g_me.d_YHYMatr, &g_me.c_YHYMatr, (size_t)nn*nn*T);
    me_grow_d(&g_me.d_YHYCoeff,&g_me.c_YHYCoeff,(size_t)T);
    me_grow_d(&g_me.d_massmat, &g_me.c_massmat, (size_t)nn*nn);
    me_grow_d(&g_me.d_charge,  &g_me.c_charge,  (size_t)nn);
    me_grow_d(&g_me.d_H,       &g_me.c_H,       (size_t)npairs);
    me_grow_d(&g_me.d_S,       &g_me.c_S,       (size_t)npairs);
    me_grow_d(&g_me.d_Dk,      &g_me.c_Dk,      (size_t)npairs*npt2);
    me_grow_d(&g_me.d_Dl,      &g_me.c_Dl,      (size_t)npairs*npt2);
    me_grow_i(&g_me.d_klist,   &g_me.c_klist,   (size_t)npairs);
    me_grow_i(&g_me.d_llist,   &g_me.c_llist,   (size_t)npairs);
    me_grow_i(&g_me.d_grfl,    &g_me.c_grfl,    (size_t)npairs);

    // d_Dl entries of non-grad_l pairs are never written by the kernel, so they
    // must be zeroed here. d_H/d_S/d_Dk are fully overwritten -> no memset.
    CUDA_CHECK(cudaMemset(g_me.d_Dl, 0, (size_t)npairs*npt2*sizeof(double)));

    // Per-call inputs (change every call): always uploaded.
    CUDA_CHECK(cudaMemcpy(g_me.d_NLP,   NonlinParam, (size_t)Nmax*nnp*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(g_me.d_klist, k_list,      (size_t)npairs*sizeof(int),      cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(g_me.d_llist, l_list,      (size_t)npairs*sizeof(int),      cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(g_me.d_grfl,  grad_l_list, (size_t)npairs*sizeof(int),      cudaMemcpyHostToDevice));

    // Run-constant inputs: uploaded once (every call in legacy mode).
    if (!g_me.constants_uploaded) {
        CUDA_CHECK(cudaMemcpy(g_me.d_YHYMatr,  YHYMatr,  (size_t)nn*nn*T*sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(g_me.d_YHYCoeff, YHYCoeff, (size_t)T*sizeof(double),       cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(g_me.d_massmat,  massmat,  (size_t)nn*nn*sizeof(double),   cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(g_me.d_charge,   charge,   (size_t)nn*sizeof(double),      cudaMemcpyHostToDevice));
        g_me.constants_uploaded = 1;
    }

    // Shared memory: H(1) + S(1) + Dk(2np) + Dl(2np) = 2 + 4*np doubles
    size_t smem = (2 + 4*nnp) * sizeof(double);
    int threads = min(T, 1024);
    matelem_batch_grad_kernel<<<npairs, threads, smem>>>(
        g_me.d_NLP, g_me.d_klist, g_me.d_llist, g_me.d_grfl,
        g_me.d_YHYMatr, g_me.d_YHYCoeff,
        T, nn, nnp, *NumOfProcs_in, *ProcID_in,
        *piraised3n2, g_me.d_massmat, g_me.d_charge, *charge0,
        *attr_scale, *rep_scale, *rep_plus, *rep_minus,
        g_me.d_H, g_me.d_S, g_me.d_Dk, g_me.d_Dl
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(Hout,  g_me.d_H,  npairs*sizeof(double),               cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Sout,  g_me.d_S,  npairs*sizeof(double),               cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Dkout, g_me.d_Dk, (size_t)npairs*npt2*sizeof(double),  cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Dlout, g_me.d_Dl, (size_t)npairs*npt2*sizeof(double),  cudaMemcpyDeviceToHost));

    if (!persist) me_release_buffers();
    if (profile) { g_prof_grad += me_now() - t0; g_calls_grad++; }
}

} // extern "C"
