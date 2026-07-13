#!/usr/bin/env python3
"""Correctness anchor for the fixed-basis GPU path (ECG_VERIFY=1).

Reads verify_dump.bin (Fortran unformatted: cbs:int32, lambda:f64, H:f64[cbs,cbs]
col-major, S:f64[cbs,cbs], x:f64[cbs]) and reports the generalized eigenpair
residual and S-normalization -- i.e. VERIFIES the eigenpair solves the system,
rather than only checking that two energies agree.

Usage:
  verify_dump.py <dump.bin>                 # residual + x^T S x + cond(S)
  verify_dump.py <cpu.bin> <gpu.bin>        # + element-wise H/S diff and dlambda
"""
import sys
import numpy as np

def load(path):
    try:
        from scipy.io import FortranFile
    except ImportError:
        sys.exit("needs scipy (pip install scipy) for Fortran-record parsing")
    f = FortranFile(path, 'r')
    cbs = int(f.read_ints(np.int32)[0])
    lam = float(f.read_reals(np.float64)[0])
    H = f.read_reals(np.float64).reshape((cbs, cbs), order='F')
    S = f.read_reals(np.float64).reshape((cbs, cbs), order='F')
    x = f.read_reals(np.float64)
    f.close()
    return cbs, lam, H, S, x

def report(tag, cbs, lam, H, S, x):
    Hx, Sx = H @ x, S @ x
    r = Hx - lam * Sx                       # generalized residual H x - lambda S x
    xSx = float(x @ Sx)                     # S-normalization (1 if x is S-normalized)
    rn = np.linalg.norm(r)
    scale = np.linalg.norm(H, 1) * np.linalg.norm(x)   # relative-residual denominator
    symH = np.max(np.abs(H - H.T))
    symS = np.max(np.abs(S - S.T))
    print(f"[{tag}] K={cbs}  lambda={lam:.12f}")
    print(f"    ||Hx - lambda Sx||_2      = {rn:.3e}")
    print(f"    relative residual         = {rn/scale:.3e}   (/(||H||_1 ||x||))")
    print(f"    x^T S x                   = {xSx:.12f}   (dev from 1: {abs(xSx-1):.2e})")
    print(f"    max|H-H^T|, max|S-S^T|    = {symH:.2e}, {symS:.2e}")
    return lam, H, S

def main():
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    lamA, HA, SA = report("A", *load(sys.argv[1]))
    if len(sys.argv) == 3:
        lamB, HB, SB = report("B", *load(sys.argv[2]))
        print("--- A vs B (e.g. CPU vs GPU matrix elements) ---")
        print(f"    |lambda_A - lambda_B|     = {abs(lamA-lamB):.3e}")
        print(f"    max|H_A - H_B|            = {np.max(np.abs(HA-HB)):.3e}")
        print(f"    max|S_A - S_B|            = {np.max(np.abs(SA-SB)):.3e}")

if __name__ == "__main__":
    main()
