# Input File Format

Table of Contents:

- [Basic structure of the input file](#1-basic-structure-of-the-input-file)
- [Input file sections](#2-input-file-sections)
  - [Header](#header)
  - [Command list](#command-list)
  - [History](#history)
  - [Basis functions](#basis-functions)
- [Creating input files](#3-creating-input-files)
- [Restarting calculations that generate a basis](#4-restarting-calculations-that-generate-a-basis)

## Basic structure of the input file

The input file `inout.txt` consists of four sections: header, list of commands, calculation history, and parameters of basis functions:

| Section | Short Description |
| :--- | :--- |
| **Header** | Defines the quantum system - number of particles, their charges, masses, permutational symmetry. Also contains the current size of the basis, current energy, parameters for the eigensolver and its statistics, as well as parameters of the distribution used in stochastic selection of variational parameters. |
| **Command List** | Contains a list of actions that need to be performed. |
| **History** | Lists the energy values that were obtained at each basis size when the basis was being built, as well as the number of optimization steps that were performed (e.g. number of objective function evaluations, number of optimization cycles). |
| **Basis functions** | Contains parameters of basis functions, one function per line. |
| | |

When calculations begin from scratch, only the header and command list need to be supplied in the input file.

Note that in any calculation that involves generation of a basis (either from scratch or by extending an existing basis) the input file also serves as an output file (hence the name - `inout.txt`). This file is updated with any incremental changes of the basis. The frequency with which the file gets overwritten can be controlled by the user.

In many practical situations one must save intermediate basis sets. This can be easily achieved with a command `SAVE_FILE`, which is described below.

## Input file sections

### Header

A typical header looks as follows (it is an example for the ground state calculations of lithium atom that begin from scratch - no basis is yet generated):

```
 BASIS_TYPE RG_0S
 PARTICLES    4
 MASSES    12786.39227775 1.0 1.0 1.0
 CHARGES   3.0 -1.0 -1.0 -1.0
 SYMMETRY  (1+P23)(1-P24)
 BASIS_SIZE  0
 CURRENT_ENERGY    1.0E30
 WHICH_EIGENVALUE  1    
 EIGVAL_TOLERANCE  1.0E-12
 INVITPARAMETER    1.000001000E0
 LAST_EIGVAL_TOL   1.0E0
 BEST_EIGVAL_TOL   1.0E0
 WORST_EIGVAL_TOL  1.0E-30
 GENERATOR_PARAM   0.7E0  1.0E0  3.0E0
```

At its end, the header includes an excerpt of the last line in the history, which may look like this:

```
 ==============================
    100 -0.7477376516187604E+01     10      0      0
 ==============================
```

This line lists the current state of the optimization process.

Note that each line begins with a keyword. Some lines in the header are optional and can be skipped. However, most are required. Lines with the keywords must appear in a specific order and no blank lines are allowed.
The description of all header keywords, in the order they appear in input file, is provided in the table below. Examples of complete input file lines for each keyword are also provided.

| Keyword | Requirement | Description |
| :--- | :---: | :--- |
| `BASIS_TYPE` | optional | Specifies the type of basis. Possible values are `CG_0S`, `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`. While this line is optional, it is highly recommended to include it in the input file to avoid confusion when running multiple calculations that use different basis types. <br> *Example* : `BASIS_TYPE RG_2D` |
| `PARTICLES` | required | Specifies the number of particles in the system. For example, for Li atom it is 4 (a nucleus + three electrons). <br> *Example* : `PARTICLES 5` |
| `FIXED_INDEX` | optional | This optional line can only be used in the case of `RG_1P` basis type. The keyword is followed by the fixed value of the z-index that should be used in calculations. <br> *Example* : `FIXED_INDEX 2` |
| `FIXED_INDEX_1` | optional | This optional line can only be used in the case of `RG_2D` or `RG_2P` basis types. The keyword is followed by the fixed value of the first integer index that should be used in calculations. <br> *Example* : `FIXED_INDEX_1 2` |
| `FIXED_INDEX_2` | optional | This optional line can only be used in the case of `RG_2D` or `RG_2P` basis types. The keyword is followed by the fixed value of the second integer index that should be used in calculations. <br> *Example* : `FIXED_INDEX_2 3` |
| `MASSES` | required | Specifies masses of all particles (in atomic units, so that mass of the electron is 1.0), starting with the first one. The values are separated by one or more spaces. For calculations of atoms or ions with a clamped nucleus the nuclear mass should be set to a very large value, e.g. 1.0E+30. <br> *Example* : `MASSES 1.0E30 1.0 1.0 1.0 1.0` |
| `CHARGES` | required | Specifies Coulomb charges of all particles (in atomic units), starting with the first one. The values are separated by one or more spaces. It is possible to use noninteger charges. <br> *Example* : `CHARGES 4.0 -1.0 -1.0 -1.0 -1.0` |
| `REPULSION_SCALING_PARAM` | optional | |
| `REPULSION_SCALING_PARAM_PLUS` | optional | |
| `REPULSION_SCALING_PARAM_MINUS` | optional | |
| `ATTRACTION_SCALING_PARAM` | optional | |
| `SYMMETRY` | required | The keyword is followed by a string that defines the Young projection operator for systems with identical particles. ECGPACK adopts spin-free formalism in the nonrelativistic calculations of the energy. The Young operator is a product of symmetrizers over all rows and antisymmetrizers over all columns of a suitable Young tableau. For the ground doublet state of lithium, a suitable Young diagram is $[ 2 \ 1]$ , which yields a symmetrizer over the first row `(1+P23)` and an antisymmetrizer over the first column `(1-P24)`. Here `Pij` denotes a pair permutation (transposition) of particles i and j. The Young operator becomes progressively more complicated for systems with a larger number of identical particles, but the structure in the form of the symmetrizer and antisymmetrizer (or, alternatively, in reverse order) remains. For example, for the doublet states of boron (a six-particle system; nucleus is particle 1, electrons are particles 2 through 6) with a Young diagram $[ 2^2 \ 1 ]$  a suitable expression for the Young operator is `(1+P23)(1+P45)(1-P24)(1-P26-P46)(1-P35)`. For the quartet states of boron with a Young diagram $[ 2 \ 1^3]$ the expression for the Young operator is `(1+P23)(1-P24)(1-P25-P45)(1-P26-P46-P56)`. If there is more than one set of identical particles in the system, the total Young operator is a product of Young operators for each set. For some systems there may be additional symmetry - e.g. charge conjugation - that can be enforced by choosing a proper expression for the Young operator. For example, for the ground state of the positronium molecule, where the first two particles are electrons and the third and fourth particles are positrons, a suitable expression is `(1+P13P24)(1+P34)(1+P12)`. <br> *Example* : `SYMMETRY (1+P23)(1-P24)(1-P25-P45)` |
| `BASIS_SIZE` | required | Indicates the current size of the basis. In new calculations it should be set to zero. <br> *Example* : `BASIS_SIZE  100` |
| `CURRENT_ENERGY` | required | The value of the current energy for the basis stored in the input file. In new calculations that start from scratch this should be set to some very large value, e.g. 1.0E30. <br> *Example* : `CURRENT_ENERGY -0.7477376516187604E+01` |
| `WHICH_EIGENVALUE` | required | Specifies which energy eigenvalue should be targeted in the calculations. 1 means the lowest, 2 means the second lowest, etc. Note that this parameter matters only if the LAPACK eigenvalue solver is used (see the next subsection for more details on that). If the inverse iteration eigenvalue solver is used this parameter can be anything. Also note that if the LAPACK eigensolver is used the value specified by `WHICH_EIGENVALUE` cannot exceed the current size of the basis. Therefore, when initiating calculations of the $k$-th excited state one must first generate or copy a basis containing at least $k+1$ ECG basis functions for the ground (or some lower excited state) and then switch to the eigenvalue number $k$ . <br> *Example* : `WHICH_EIGENVALUE 2` |
| `EIGVAL_TOLERANCE` | required | Specifies the desired accuracy of the inverse iteration eigenvalue solver. It is recommended to make `EIGVAL_TOLERANCE` as small as possible in each particular calculation but not too close to the value of the machine epsilon. In practice, due to rounding errors and some degree of linear dependencies in the generated bases, the value of `EIGVAL_TOLERANCE` should be a couple of orders of magnitude larger than the machine epsilon. Otherwise the inverse iteration eigenvalue solver may fail too often. For calculations with double precision a reasonable starting choice for `EIGVAL_TOLERANCE` is somewhere around 1.0E-12 and as the basis generation progresses it may need to be adjusted to balance between accuracy and stability. When the LAPACK eigenvalue solver is invoked in calculations, the value `EIGVAL_TOLERANCE` is not referenced and does not affect anything. <br> *Example* : `EIGVAL_TOLERANCE 1.0E-12` |
| `INVITPARAMETER` | required | |
| `LAST_EIGVAL_TOL` | required | |
| `BEST_EIGVAL_TOL` | required | |
| `WORST_EIGVAL_TOL` | required | |
| `GENERATOR_PARAM` | required | |
| | | |

### Command list

### History

### Basis functions

## Creating input files

When a large basis needs to be generated it requires multiple actions (e.g. basis enlargements `BASIS_ENL` and cyclic optimizations `OPT_CYCLE`) that need to be repeated routinely. It may be tedious to manually write the sequence of corresponding commands and their multiple arguments. For this reason it is convenient to create a simple utility that creates a desired list of commands, or, even better, prepares a ready-to-use input file `inout.txt` to initiate calculations from scratch.

Directory `ecgpack/utilities/input_file_manipulation` contains a Python script `ecg_input_file_generator.py` that can be invoked to generate an input file for a desired system with a long list of repeated actions. One can modify this script as needed.

## Restarting calculations that generate a basis
