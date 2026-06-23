# Compilation and execution

## Compilation

To compile any selected code (e.g. `RG_0S`) using a specific toolchain and options one can go to the code directory (e.g. `RG_0S/`) and run the `make` command with the corresponding arguments. Please look at the header of the `Makefile` in the code directory as it provides a brief explanation of the arguments.

### Precision
Note each code in this collection can be compiled with a different working precision specified by the user (provided there is compiler and hardware architecture support) as an argument in the `make` command: 

* `PREC=8` - fp64 (double precision)
* `PREC=10` - fp80 (extended precision), hardware support available on x86 CPUs only
* `PREC=16` - fp128 (quadruple precision), emulation mode

Because calculations with fp128 floating-point numbers use emulation on pretty much any modern CPU, they are very slow. As a rule of thumb, the speed of computations and memory allocation requirements of the codes in this collection on x86 CPUs scale as follows:

| Precision | CPU time (relative) | Memory amount (relative) |
| --- | --- | --- |
| fp64  | 1 | 1 |
| fp80  | 5 | 2 |
| fp128 | 100 | 2 |
| | | |

Therefore, the quadruple precision (roughly a factor of ~100 slower than double precision) should be used only in special cases where it is absolutely necessary and feasible.  

### Linear algebra library

The energy codes can be linked against different BLAS/LAPACK implementations, selected with the `LINALG` argument of the `make` command:

* `LINALG=netlib` (default) - the bundled, lightly modified netlib reference BLAS/LAPACK is compiled from source (`src/BLAS.f`/`src/LAPACK.f`); no external library is required
* `LINALG=mkl` - Intel Math Kernel Library (MKL)
* `LINALG=lblas` - an optimized BLAS/LAPACK exposed through the `-lblas`/`-llapack` symbolic links
* `LINALG=openblas` - OpenBLAS
* `LINALG=aocl` - AMD Optimizing CPU Libraries (AOCL-BLAS and AOCL-LAPACK)

An optimized library can be selected only for `PREC=8`; for `PREC=10` and `PREC=16` only `LINALG=netlib` is available, since optimized BLAS/LAPACK libraries do not support extended or quadruple precision. In the off-diagonal matrix-element codes the `LINALG` argument is accepted but has no effect, as those codes do not link an external BLAS/LAPACK library.

### Number of particles

Note that the maximum number of particles in the calculations is compiled in; it is not a runtime argument. A code compiled for a certain number of particles will not execute with input files where the number of particles is larger. For best performance, it is **strongly recommended** that for any calculation of a system with $N$ particles one uses a corresponding binary compiled for exactly the same maximum number of particles.

### Batch compilation of multiple code variants 

The easiest way to compile all or some number of selected codes in **one step** on a specific machine/OS using specific toolchains, precision, etc. is to invoke the `build.bash` script located in the root directory. This script requires arguments. **Please read its source or run it with no arguments** to see instructions regarding how to run it properly. When this script is run, it will save all individual binaries in directory `ecgpack/bin/<toolchainname>/<configuration>`, where `<toolchainname>` (e.g. `systemdefault`) is the name of the toolchain specified and `<configuration>` can be either `debug` (slow and unoptimized binary suitable for debugging) or `release` (fast and optimized binary suitable for production work). The binary files will be named `<code_name>_N<particle_number>_P<precision>_<linalg>` (e.g. `RG_0S_N4_P8_netlib` - `RG_0S` code, 4 particles, double precision, bundled netlib BLAS/LAPACK). The `<linalg>` suffix records the selected linear algebra library (`netlib`, `mkl`, `lblas`, `openblas`, or `aocl`) and is always present.

## Execution

Any compiled binary file can be executed using multiple MPI processes in a standard way:

```bash
mpirun -np <NPROCS> <BINARYFILE>
```

For example, suppose a user has compiled all codes with a help of `build.bash` script (invoking the `systemdefault` toolchain, which on most Linux machines means `gfortran` compiler with `OpenMPI`). Then the user creates a work directory `ecgpack/jobs/test` for a test calculation of a 4-particle system with the `RG_0S` ($L=0$) basis and places a relevant input file `inout.txt` in that directory. To launch a test calculation with double precision using 12 CPU cores the execution command should be
```bash
mpirun -np 12 ../../bin/systemdefault/release/RG_0S_N4_P8_netlib
```

Note that the input file `inout.txt` (as well as other relevant files, if any) should be located in the same directory where the execution of the code takes place.