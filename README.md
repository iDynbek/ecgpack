# NOTES ON BUILDING AND USING THE ECG CODE

## Short description
ECG is a collection of closely related computer codes for performing variational calculations of few-electron / few-body quantum systems using all-particle explicitly correlated Gaussian (ECG) basis functions that have been developed in the group of Sergiy Bubin (Physics Department, Nazarbayev University).

## Notations
Parts of the description below use formulas in the LaTeX format (sandwiched between dollar signs). For comfortable reading it is advised to use a Markdown viewer that can render these formulas in a human-readable form, such as Microsoft Visual Studio Code. The conventions for the mathematical notations are adopted from: 
* T. Shomenov and S. Bubin, Phys. Rev. E 108, 065308 (2023) [link](https://doi.org/10.1103/PhysRevE.108.065308)
* S. Bubin and L. Adamowicz, J. Chem. Phys. 128, 114107 (2008) [link](https://doi.org/10.1063/1.2894866)
* S. Bubin and L. Adamowicz, J. Chem. Phys. 124, 224317 (2006) [link](https://doi.org/10.1063/1.2204605)

## Directory structure
The following folders exist within the repository. The first two letters in the name (`RG` or 'CG`) stand for the real or complex Gaussians, the other letters relate to the orbital angular momenta of individual electrons and the total orbital angular momentum $L$.

* `archive` - An archive of some older versions of the ECG codes. This should normally not be used by anyone, unless you know what you are doing.

* `CG_0S` - The code for energy and wavefunction calculations with complex spherically symmetric gaussians ($L=0$, even parity) that have the form $$\phi_k=\exp [ \bm{r}' (C_k \otimes \mathbf{I}) \bm{r}]=\exp [ \bm{r}' ((A_k+i B_k) \otimes \mathbf{I}) \bm{r}].$$
The CG_0S code is currently work in progress and lags behind the corresponding RG_0S code in terms of features and the quality of the implementation.  

* `RG_0S` - The code for energy and wavefunction calculations with real spherically symmetric gaussians ($L=0$, even parity) that have the form $$\phi_k=\exp [ \bm{r}' (A_k \otimes \mathbf{I}) \bm{r}].$$

* `RG_1P` - The code for energy and wavefunction calculations with real ECGs ($L=1$, odd parity) that have the form $$\phi_k=z_{i_k}\exp [ \bm{r}' (A_k \otimes \mathbf{I}) \bm{r}].$$

* `RG_2D` - The code for energy and wavefunction calculations with real ECGs ($L=2$, even parity) that have the form $$\phi_k=(x_{i_k} x_{j_k} + y_{i_k} y_{j_k} - 2 z_{i_k} z_{j_k} ) \exp [ \bm{r}' (A_k \otimes \mathbf{I}) \bm{r}],$$ where the integer index $i_k$ can be either different or the same as $j_k$.

* `RG_2P` - The code for energy and wavefunction calculations with real ECGs ($L=1$, even parity) that have the form $$\phi_k=(x_{i_k} y_{j_k} - x_{j_k} y_{i_k}) \exp [ \bm{r}' (A_k \otimes \mathbf{I}) \bm{r}].$$

* `RG_0S-1P` - The code for the calculation of the offdiagonal matrix elements between states described with RG_0S and RG_1P bases. Currently this includes the evaluation of the transition dipole moments in the length and velocity gauges.

* `RG_0S-2P` - The code for the calculation of the offdiagonal matrix elements between states described with RG_0S and RG_2P bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions.

* `RG_0S-2D` - The code for the calculation of the offdiagonal matrix elements between states described with RG_0S and RG_2D bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions.

* `RG_1P-2D` - The code for the calculation of the offdiagonal matrix elements between states described with RG_1P and RG_2D bases. Currently this includes the evaluation of the transition dipole moments in the length and velocity gauges.

* `RG_1P-1P` - The code for the calculation of the offdiagonal matrix elements between states described with RG_1P and RGL_1P bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions.

* `RG_2D-2D` - The code for the calculation of the offdiagonal matrix elements between states described with RG_2D and RG_2D bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions.

* `RG_2P-2D` - The code for the calculation of the offdiagonal matrix elements between states described with RG_2P and RG_2D bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions.

* `RG_2P-2P` - The code for the calculation of the offdiagonal matrix elements between states described with RG_2P and RG_2P bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions.

Each folder with a code, namely `CG_0S`, `RG_S0`, `RG_1P`, `RG_2D`, `RG_2P`, `RG_0S-1P`, `RG_0S-2P`, `RG_1P-2D`, `RG_2D-2D`,`RG_2P-2D`, and `RG_2P-2P` has a `Makefile` and a subfolder `src` that contains the actual source. The source structure and the structure of the makefiles are very similar for all codes. In fact, some of the makefiles are identical.

## Compilation
The easiest way to compile all or some selected codes on a specific machine/OS using specific toolchains, precision, etc. is to invoke the `build.bash` script located in the repository's root directory. This script requires arguments. Please read its source or run it with no arguments to see instructions regarding how to use it properly.

## Execution

## Description of script commands
