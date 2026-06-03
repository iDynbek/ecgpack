# NOTES ON COMPILING AND USING THE ECG CODES

## Short description

ECG is a collection of closely related computer codes for performing variational calculations of quantum few-body systems (few-electron atoms and molecules) using all-particle explicitly correlated Gaussian (ECG) basis functions that have been developed by the research group of Sergiy Bubin (Physics Department, Nazarbayev University) and collaborators.

## Notations

Mathematical formulas in these notes are written in the LaTeX format (sandwiched between dollar signs). For comfortable reading it is advised to use a Markdown viewer that can render these formulas in a human-readable form, such as Microsoft Visual Studio Code. The conventions for the mathematical notations are adopted from:

* [T. Shomenov and S. Bubin, Phys. Rev. E 108, 065308 (2023)](https://doi.org/10.1103/PhysRevE.108.065308)
* [S. Bubin and L. Adamowicz, J. Chem. Phys. 128, 114107 (2008)](https://doi.org/10.1063/1.2894866)
* [S. Bubin and L. Adamowicz, J. Chem. Phys. 124, 224317 (2006)](https://doi.org/10.1063/1.2204605)

## Directory structure

The ECG project repository has the following directory structure:

| Directory Tree | Description |
| -------- | -------- |
| `ecg/` | Root directory |
| `├──[archive]/` | An archive of some older versions of the ECG codes. This should normally not be used or looked at by anyone, unless you know what you are doing. |
| `├──[bin]/` | Binary files for calculations (created by user) |
| `├──[CG_OS]/` | The code for energy and wavefunction calculations with complex spherically symmetric gaussians ($L=0$, even parity) that have the form $$\phi_k=\exp [ \mathbf{r}' (C_k \otimes \mathbf{I}) \mathbf{r}]=\exp [ \mathbf{r}' ((A_k+i B_k) \otimes \mathbf{I}) \mathbf{r}].$$ The CG_0S code is currently work in progress and lags behind the corresponding RG_0S code in terms of features and the quality of the implementation. |
| `├──[jobs]/` | Work directory for calculations (created by user) |
| `├──[RG_OS]/` | The code for energy and wavefunction calculations with real spherically symmetric gaussians ($L=0$, even parity) that have the form $$\phi_k=\exp [ \mathbf{r}' (A_k \otimes \mathbf{I}) \mathbf{r}].$$ |
| `├──[RG_OS-1P]/` | The code for the calculation of the transition electric dipole moments for states that are expanded using RG_0S and RG_1P bases. Currently this includes the transition dipole moment operators in the length and velocity gauges. |
| `├──[RG_OS-2D]/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_0S and RG_2D bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `├──[RG_OS-2P]/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_0S and RG_2P bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `├──[RG_1P]/` | The code for energy and wavefunction calculations with real ECGs ($L=1$, odd parity) that have the form $$\phi_k=z_{i_k}\exp [ \mathbf{r}' (A_k \otimes \mathbf{I}) \mathbf{r}].$$ |
| `├──[RG_1P-1P]/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_1P and RGL_1P bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `├──[RG_1P-2D]/` | The code for the calculation of the transition electric dipole moments for states that are expanded using RG_1P and RG_2D bases. Currently this includes the transition dipole moment operators in the length and velocity gauges. |
| `├──[RG_1P-2P]/` | The code for the calculation of the transition electric dipole moments for states that are expanded using RG_1P and RG_2P bases. Currently this includes the transition dipole moment operators in the length and velocity gauges. |
| `├──[RG_2D]/` | The code for energy and wavefunction calculations with real ECGs ($L=2$, even parity) that have the form $$\phi_k=(x_{i_k} x_{j_k} + y_{i_k} y_{j_k} - 2 z_{i_k} z_{j_k} ) \exp [ \mathbf{r}' (A_k \otimes \mathbf{I}) \mathbf{r}],$$ where the integer index $i_k$ can be either different or the same as $j_k$. |
| `├──[RG_2D-2D]/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_2D and RG_2D bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `├──[RG_2P]/` | The code for energy and wavefunction calculations with real ECGs ($L=1$, even parity) that have the form $$\phi_k=(x_{i_k} y_{j_k} - x_{j_k} y_{i_k}) \exp [ \mathbf{r}' (A_k \otimes \mathbf{I}) \mathbf{r}].$$ |
| `├──[RG_2P-2D]/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_2P and RG_2D bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `├──[RG_2P-2P]/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_2P and RG_2P bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `└──[utilities]/` | Various utilities |

Each directory with a code contains a `Makefile` and a subdirectory `src` with the actual source. The source structure and the structure of the makefiles are very similar for all codes. In fact, some of the makefiles are identical.

## Compilation

The easiest way to compile all or some selected codes on a specific machine/OS using specific toolchains, precision, etc. is to invoke the `build.bash` script located in the root directory. This script requires arguments. Please read its source or run it with no arguments to see instructions regarding how to use it properly.

## Execution

## Input file
