# $S^e \rightarrow P^o$ Transition Dipole Moment Calculation

## Description

This calculation evaluates the **electric-dipole transition matrix element** between an initial $S^e$ state and a final $P^o$ state using the ECG transition code.

The ECG code calculates only the **transition dipole moment**. It does **not** calculate the line strength or oscillator strength directly.

The line strength and oscillator strength are calculated afterward in a separate post-processing step, for example in an Excel worksheet, using the calculated transition dipole moment, transition energy, and the required angular-momentum factors.

Only the wave-function files of the initial and final states are required when the transition code is executed. No input file such as `inout.txt` is required at that stage.

The initial and final wave functions must correspond to the **same isotope** and must use the same particle masses and particle ordering. Otherwise, the transition dipole moment calculation will fail or the resulting matrix element will not be physically meaningful.

## 1. Prepare the initial and final wave functions

Two independently optimized ECG wave functions are required:

| State | Role in the transition |
|---|---|
| $S^e$ state | Initial state |
| $P^o$ state | Final state |

The wave functions must be calculated for the same isotope and with consistent physical parameters.

The calculation files are organized in the following three folders:

```text
intial_state_Li-2Se/
final_state_Li-2Po/
transition_dipole_moment/
```

The folder contents are organized as follows:

- `intial_state_Li-2Se` contains the input files for the initial $S^e$-state calculation and the generated initial-state wave function.
- `final_state_Li-2Po` contains the input files for the final $P^o$-state calculation and the generated final-state wave function.
- `transition_dipole_moment` contains the wave-function files required by the transition code and the resultant output file.

## 2. Save each wave function

Add the following instruction to the corresponding `inout.txt` file used in each ECG state calculation:

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

Copy the generated wave functions from the initial- and final-state folders to the `transition_dipole_moment` folder and rename them accordingly.

For example, when the three folders are located in the same parent directory:

```bash
cd transition_dipole_moment
cp ../intial_state_Li-^2S^e/wavefunction.txt wf_state0.txt
cp ../final_state_Li-^2P^o/wavefunction.txt wf_state1.txt
```

## 4. Prepare the working directory

The `transition_dipole_moment` folder is the working directory for the transition calculation. Place the following files in this folder:

```text
RG_0S-1P
wf_state0.txt
wf_state1.txt
```

An `inout.txt` file is not required in this directory because the transition code reads the two wave-function files directly.

Before the calculation, the directory should therefore contain at least:

```text
transition_dipole_moment/
├── RG_0S-1P
├── wf_state0.txt
└── wf_state1.txt
```

After the calculation, the resultant output file is also generated in the `transition_dipole_moment` folder.

## 5. Run the transition calculation

Run the `RG_0S-1P` executable from the `transition_dipole_moment` working directory.

For example:

```bash
./RG_0S-1P
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

This conversion, together with the subsequent calculation of the line strength and oscillator strength, is performed separately from the ECG transition calculation.

For example, these calculations may be carried out in an Excel worksheet using:

- the transition dipole moment obtained from `RG_0S-1P`;
- the transition energy;
- the statistical weight of the initial fine-structure level;
- the required $3j$-symbol factor;
- the required $6j$-symbol factor; and
- the length- or velocity-gauge oscillator-strength expression.

## Calculation workflow

```text
Calculate the initial S-state wave function in
intial_state_Li-^2S^e
                         +
Calculate the final P-state wave function in
final_state_Li-^2P^o
                         ↓
Save both wave functions using SAVE_HSWF
                         ↓
Copy and rename the initial S-state wave function as
transition_dipole_moment/wf_state0.txt
                         ↓
Copy and rename the final P-state wave function as
transition_dipole_moment/wf_state1.txt
                         ↓
Run RG_0S-1P in the transition_dipole_moment folder
                         ↓
Generate the resultant output file and obtain
the transition dipole moment
                         ↓
Calculate the line strength and oscillator strength separately,
for example in Excel
```

## Important notes

1. `RG_0S-1P` calculates the transition dipole moment, not the oscillator strength.
2. The initial state must be the $S^e$ state stored in `wf_state0.txt`.
3. The final state must be the $P^o$ state stored in `wf_state1.txt`.
4. Both wave functions must correspond to the same isotope.
5. Both wave functions must use the same particle ordering and masses.
6. No `inout.txt` file is needed when `RG_0S-1P` is executed.
7. The oscillator strength must be calculated afterward in a separate post-processing step.
