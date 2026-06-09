# NOTES ON COMPILING AND USING THE ECG CODES

## Short description

ECG is a collection of closely related parallel computer codes for performing variational calculations of quantum few-body systems (few-electron atoms and molecules) using all-particle explicitly correlated Gaussian (ECG) basis functions that have been developed by the research group of Sergiy Bubin (Physics Department, Nazarbayev University) and collaborators. The codes are written in Fortran and use MPI for parallelism.

## Notations

Mathematical formulas in these and other notes are written in the LaTeX format (sandwiched between dollar signs). For comfortable reading it is advised to use a Markdown viewer that can render these formulas in a human-readable form, such as Microsoft Visual Studio Code. Basic theory as well as conventions for the mathematical notations are available in the following references:

* [T. Shomenov and S. Bubin, Phys. Rev. E 108, 065308 (2023)](https://doi.org/10.1103/PhysRevE.108.065308)
* [S. Bubin and L. Adamowicz, J. Chem. Phys. 128, 114107 (2008)](https://doi.org/10.1063/1.2894866)
* [S. Bubin and L. Adamowicz, J. Chem. Phys. 124, 224317 (2006)](https://doi.org/10.1063/1.2204605)

## Directory structure

The ECG project repository has the following directory structure:

| Directory | Description |
| -------- | -------- |
| `ecg/` | Root directory |
| `ecg/archive/` | An archive of some older versions of the ECG codes. This should normally not be used or looked at by anyone, unless you know what you are doing. |
| `ecg/bin/` | Binary files for calculations (created by user) |
| `ecg/CG_0S/` | The code for energy and wavefunction calculations with complex spherically symmetric Gaussians ($L=0$, even parity) that have the form $$\phi_k=\exp [ \mathbf{r}' (C_k \otimes \mathbf{I}) \mathbf{r}]=\exp [ \mathbf{r}' ((A_k+i B_k) \otimes \mathbf{I}) \mathbf{r}].$$ The CG_0S code is currently a work in progress and lags behind the corresponding RG_0S code in terms of features and the quality of the implementation. |
| `ecg/jobs/` | Work directory for calculations (created by user) |
| `ecg/RG_0S/` | The code for energy and wavefunction calculations with real spherically symmetric Gaussians ($L=0$, even parity) that have the form $$\phi_k=\exp [ \mathbf{r}' (A_k \otimes \mathbf{I}) \mathbf{r}].$$ |
| `ecg/RG_0S-1P/` | The code for the calculation of the transition electric dipole moments for states that are expanded using RG_0S and RG_1P bases. Currently this includes the transition dipole moment operators in the length and velocity gauges. |
| `ecg/RG_0S-2D/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_0S and RG_2D bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `ecg/RG_0S-2P/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_0S and RG_2P bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `ecg/RG_1P/` | The code for energy and wavefunction calculations with real ECGs ($L=1$, odd parity) that have the form $$\phi_k=z_{i_k}\exp [ \mathbf{r}' (A_k \otimes \mathbf{I}) \mathbf{r}].$$ |
| `ecg/RG_1P-1P/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_1P and RGL_1P bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `ecg/RG_1P-2D/` | The code for the calculation of the transition electric dipole moments for states that are expanded using RG_1P and RG_2D bases. Currently this includes the transition dipole moment operators in the length and velocity gauges. |
| `ecg/RG_1P-2P/` | The code for the calculation of the transition electric dipole moments for states that are expanded using RG_1P and RG_2P bases. Currently this includes the transition dipole moment operators in the length and velocity gauges. |
| `ecg/RG_2D/` | The code for energy and wavefunction calculations with real ECGs ($L=2$, even parity) that have the form $$\phi_k=(x_{i_k} x_{j_k} + y_{i_k} y_{j_k} - 2 z_{i_k} z_{j_k} ) \exp [ \mathbf{r}' (A_k \otimes \mathbf{I}) \mathbf{r}],$$ where the integer index $i_k$ can be either different or the same as $j_k$. |
| `ecg/RG_2D-2D/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_2D and RG_2D bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `ecg/RG_2P/` | The code for energy and wavefunction calculations with real ECGs ($L=1$, even parity) that have the form $$\phi_k=(x_{i_k} y_{j_k} - x_{j_k} y_{i_k}) \exp [ \mathbf{r}' (A_k \otimes \mathbf{I}) \mathbf{r}].$$ |
| `ecg/RG_2P-2D/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_2P and RG_2D bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `ecg/RG_2P-2P/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_2P and RG_2P bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `ecg/utilities/` | Various utilities |

Each directory with a code contains a `Makefile` and a subdirectory `src` with the actual source. The source structure and the structure of the makefiles are very similar for all codes. In fact, some of the makefiles for different codes in the collection are identical. Each directory with a code also contains an identical subdirectory `.vscode` with standard JSON configuration files for Microsoft Visual Studio Code (VSCode). These are user-independent and provide the capability to build and debug each code in its directory using VSCode.

Each directory with a code also contains a subdirectory called `sample_input`. These contain simple examples of input files (files `inout.txt`) as well as instructions (in files `README.md`) on how to execute those sample calculations and any additional files or scripts.

## Root directory files

There are several files located in the project's root directory `ecg/`. Their description is provided below:

| File | Description |
| -------- | -------- |
| `README.md` | This repository manual file |
| `CLAUDE.md` | Configuration file for Anthropic Claude Code to establish persistent project context | 
| `build.bash` | A Bash-script for batch compilation of multiple codes corresponding to a different number of particles, toolchains, configurations, precision, etc. It is convenient for building a large number of different binaries that are later used in production calculations. The generated binaries are automatically moved to directory `/ecg/bin`. For more information run this script in a terminal without arguments or read its header. |
| `.code-workspace` | A JSON configuration file for Microsoft Visual Studio Code (VSCode) that contains information used to group separate code project directories into a single, unified workspace that can be opened in VSCode. For tips on using VSCode see [Use of Microsoft Visual Studio Code](#use-of-microsoft-visual-studio-code) below. |

## Compilation

To compile any selected code (e.g. `RG_0S`) using a specific toolchain and options one can go to the code directory (e.g. `RG_0S/`) and run the `make` command with the corresponding arguments. Please look at the header of the `Makefile` in the code directory as it provides a brief explanation of the arguments.

### Precision
Note each code in this collection can be compiled with a different working precision specified by the user (provided there is compiler and hardware architecture support) as an argument in the `make` command: 

* `PREC=8` - fp64 (double precision)
* `PREC=10` - fp80 (extended precision), hardware support available on x86 CPUs only
* `PREC=16` - fp128 (quadruple precision), emulation mode

Because calculations with fp128 floating-point numbers use emulation on pretty much any modern CPU, they are very slow. As a rule of thumb, the speed of computations and memory allocation requirements of the codes in this collection on x86 CPUs scale as follows:

| Precision | CPU time (relative) | Memory amount (relative)
| --- | --- | --- |
| fp64  | 1 | 1 |
| fp80  | 5 | 2 |
| fp128 | 100 | 2 |

Therefore, the quadruple precision (roughly a factor of ~100 slower than double precision) should be used only in special cases where it is absolutely necessary and feasible.  

### Number of particles

Note that the maximum number of particles in the calculations is compiled in; it is not a runtime argument. A code compiled for a certain number of particles will not execute with input files where the number of particles is larger. For best performance, it is **strongly recommended** that for any calculation of a system with $N$ particles one uses a corresponding binary compiled for exactly the same maximum number of particles.

### Batch compilation of multiple code variants 

The easiest way to compile all or some number of selected codes in **one step** on a specific machine/OS using specific toolchains, precision, etc. is to invoke the `build.bash` script located in the root directory. This script requires arguments. **Please read its source or run it with no arguments** to see instructions regarding how to run it properly. When this script is run, it will save all individual binaries in directory `ecg/bin/<toolchainname>/<configuration>`, where `<toolchainname>` (e.g. `systemdefault`) is the name of the toolchain specified and `<configuration>` can be either `debug` (slow and unoptimized binary suitable for debugging) or `release` (fast and optimized binary suitable for production work). The binary files will be named `<code_name>_N<particle_number>_P${precision}` (e.g. `RG_0S_N4_P8` - `RG_0S` code, 4 particles, double precision).

## Execution

Any compiled binary file can be executed using multiple MPI processes in a standard way:

```bash
mpirun -np <NPROCS> <BINARYFILE>
```

For example, suppose a user has compiled all codes with a help of `build.bash` script (invoking the `systemdefault` toolchain, which on most Linux machines means `gfortran` compiler with `OpenMPI`). Then the user creates a work directory `ecg/jobs/test` for a test calculation of a 4-particle system with the `RG_0S` ($L=0$) basis and places a relevant input file `inout.txt` in that directory. To launch a test calculation with double precision using 12 CPU cores the execution command should be
```bash
mpirun -np 12 ../../bin/systemdefault/release/RG_0S_N4_P8
```

Note that the input file `inout.txt` (as well as other relevant files, if any) should be located in the same directory where the execution of the code takes place.

## Input file

## Use of Microsoft Visual Studio Code
