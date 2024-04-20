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

* `RGL0` - The code for energy and wavefunction calculations with real spherically symmetric gaussians ($L=0$, even parity) that have the form $$\phi_k=\exp [ \bm{r}' (A_k \otimes \mathbf{I}) \bm{r}].$$

* `RGL1` - The code for energy and wavefunction calculations with real ECGs ($L=1$, odd parity) that have the form $$\phi_k=z_{i_k}\exp [ \bm{r}' (A_k \otimes \mathbf{I}) \bm{r}].$$

* `RGL2P` - The code for energy and wavefunction calculations with real ECGs ($L=1$, even parity) that have the form $$\phi_k=(x_{i_k} y_{j_k} - x_{j_k} y_{i_k}) \exp [ \bm{r}' (A_k \otimes \mathbf{I}) \bm{r}].$$

* `RGL2D` - The code for energy and wavefunction calculations with real ECGs ($L=2$, even parity) that have the form $$\phi_k=(x_{i_k} x_{j_k} + y_{i_k} y_{j_k} - 2 z_{i_k} z_{j_k} ) \exp [ \bm{r}' (A_k \otimes \mathbf{I}) \bm{r}],$$ where the integer index $i_k$ can be either different or the same as $j_k$.

* `CGL0` - The code for energy and wavefunction calculations with complex spherically symmetric gaussians ($L=0$, even parity) that have the form $$\phi_k=\exp [ \bm{r}' (C_k \otimes \mathbf{I}) \bm{r}]=\exp [ \bm{r}' ((A_k+i B_k) \otimes \mathbf{I}) \bm{r}].$$
The CGL0 code is currently work in progress and lags behind the corresponding RGL0 code in terms of features and the quality of the implementation.  

* `RGL01` - The code for the calculation of the offdiagonal matrix elements between states described with RGL0 and RGL1 bases. Currently this includes the evaluation of the transition dipole moments in the length and velocity gauges.

* `RGL02P` - The code for the calculation of the offdiagonal matrix elements between states described with RGL0 and RGL2P bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions.

* `RGL2P2D` - The code for the calculation of the offdiagonal matrix elements between states described with RGL2P and RGL2D bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions.
* `RGL2P2P` - The code for the calculation of the offdiagonal matrix elements between states described with RGL2P and RGL2P bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions.
* `archive` - An archive of some older versions of the ECG codes. This should normally not be used by anyone, unless you know what you are doing.

Each folder with a code, namely `RGL0`, `RGL1`, `RGL2P`, `RGL2D`, `CGL0`, `RGL01`, `RGL02P`, `RGL2P2D`, and `RGL2P2P` has a `Makefile` and a subfolder `src` that contains the actual source. The source structure and the structure of the makefiles are very similar for all codes. In fact, some of the makefiles are identical.

## Compilation

## Execution

## Description of script commands
