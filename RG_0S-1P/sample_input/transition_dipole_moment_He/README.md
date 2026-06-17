# $S^e \rightarrow P^o$ Transition Dipole Moment Calculation (Helium atom)

## Description

This calculation evaluates the **electric-dipole transition matrix element** between the initial $1~^1S^e$ state and the final $2~^1P^o$ state of the helium atom using the ECG transition dipole moment code.

The RG_0S-1P code calculates only the **transition dipole moment**. It does **not** calculate the line strength or oscillator strength directly.

The line strength and oscillator strength are calculated afterward in a separate post-processing step, for example in an Excel worksheet, using the calculated transition dipole moment, transition energy, and the required angular-momentum factors. For the definitions of the reduced transition matrix elements, line strengths, and length- and velocity-form oscillator strengths, see:

S. Nasiri, L. Adamowicz, and S. Bubin, *Physical Review A* **112**, 062809 (2025), https://doi.org/10.1103/qrtf-56np.

Only the wave-function files of the initial ($S^e$) and final ($P^o$) states are required when the transition code is executed. 
No input file such as `inout.txt` is required at that stage.

The initial and final wave functions must correspond to the **same isotope** and must use the same particle masses and particle ordering.
Otherwise, the transition dipole moment calculation will fail or the resulting matrix element will not be physically meaningful.


## 1. Required files 

Two independently optimized ECG wave functions are required:

| State | Role in the transition |
|---|---|
| $S^e$ state | Initial state |
| $P^o$ state | Final state |

The wave functions must be calculated for the same isotope and with consistent physical parameters.

The example files are organized in the following three folders:

```text
transition_dipole_moment_He
   initial_state_He-1Se/
   final_state_He-2Po/
   dipole_moment/
```

The folder contents are organized as follows:

- `initial_state_He-1Se` contains the input file for generating initial-state wave function ($1~^1S^e$).
- `final_state_He-2Po` contains the input file for generating the final-state wave function ($2~^1P^o$).
- `dipole_moment` contains the wave-function files required by the transition code.


## 2. Save each wave function

Add the following instruction to the corresponding `inout.txt` files used in each ECG state calculation:

```text
SAVE_HSWF I 100 none none none wavefunction.txt
```

The relevant arguments are:

| Argument | Description |
|---|---|
| `SAVE_HSWF` | Keyword for writing the wave function |
| `I` | Eigensolver type |
| `100` | Number of ECG basis functions |
| `wavefunction.txt` | Name of the saved wave-function file |

The basis-set size specified in the `SAVE_HSWF` instruction must be consistent with the wave function being saved.
Run the $S^e$-state and $P^o$-state calculations separately and save both wave functions.


## 3. Rename the wave-function files

The `RG_0S-1P` transition code expects the following file names:

| State | Required file name |
|---|---|
| Initial $S^e$ state | `wf_state0.txt` |
| Final $P^o$ state | `wf_state1.txt` |

Rename or copy the saved wave functions accordingly.

Copy the generated wave functions from the initial- and final-state folders to the `dipole_moment` folder and rename them accordingly.

For example, when the three folders are located in the same parent directory:

```bash
cd transition_dipole_moment_He
cd dipole_moment
cp ../initial_state_He-1Se/wavefunction.txt wf_state0.txt
cp ../final_state_He-2Po/wavefunction.txt wf_state1.txt
```


## 4. Prepare the working directory

The `dipole_moment` folder is the working directory for the transition calculation.
Place the following files in this folder:

```text
RG_0S-1P
wf_state0.txt
wf_state1.txt
```

An `inout.txt` file is not required in this directory because the transition code reads the two wave-function files directly.

Before the calculation, the directory should therefore contain at least:

```text
dipole_moment/
├── RG_0S-1P
├── wf_state0.txt
└── wf_state1.txt
```

After the calculation, the resultant output file is also generated in the `dipole_moment` folder.


## 5. Run the transition calculation

Run the `RG_0S-1P` executable from the `dipole_moment` working directory.

For example:

```bash
./RG_0S-1P >  out.txt
```

The program reads:

```text
wf_state0.txt
wf_state1.txt
```

and calculates the transition dipole moment for the $S^e\rightarrow P^o$ transition.


## 6. Use the calculated transition dipole moment

The result obtained from `RG_0S-1P` is the transition dipole moment associated with the selected component of the electric-dipole operator.

For the standard $S^e\rightarrow P^o$ calculation using the $z$ component, the directly calculated matrix element is converted to the orbital reduced transition dipole moment using the appropriate angular factor.

The orbital reduced transition matrix element is first obtained from the calculated component using the appropriate angular-momentum factor.
The line strength is then calculated from the squared magnitude of the reduced matrix element. Finally, the oscillator strength is calculated using the line strength, transition energy, and initial-state statistical weight.

For example, these calculations may be carried out in an Excel worksheet using:

- the transition dipole moment obtained from `RG_0S-1P` code;
- the transition energy;
- the statistical weight of the initial fine-structure level;
- the required $3j$-symbol factor;
- the required $6j$-symbol factor; and
- the length- or velocity-gauge oscillator-strength expression.


## Calculation workflow

```text

Save the initial S-state wave function in initial_state_He-1Se                             
                         +
Save the final P-state wave function in final_state_He-2Po
                         ↓
Copy and rename the initial S-state wave function as
dipole_moment/wf_state0.txt
                         ↓
Copy and rename the final P-state wave function as
dipole_moment/wf_state1.txt
                         ↓
Run RG_0S-1P in the dipole_moment folder
                         ↓
Generate the resultant output file and obtain
the transition dipole moment
                         ↓
Calculate the line strength and oscillator strength separately,
for example in Excel
```
