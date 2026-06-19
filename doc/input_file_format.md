# Input File Format

## Basic structure of the input file

The input file `inout.txt` consists of four sections: header, command list, history, and basis functions:

| Section | Short Description |
| --- | --- |
| Header | Defines the quantum system - number of particles, their charges, masses, permutational symmetry. Also contains the current size of the basis, current energy, parameters for the eigensolver and its statistics, as well as parameters of the distribution used in stochastic selection of variational parameters. |
| Command List | Contains a list of actions that need to be performed. |
| History | Lists the energy values that were obtained at each basis size when the basis was being built, as well as the number of optimization steps that were performed (e.g. number of objective function evaluations, number of optimization cycles). |
| Basis functions | Contains parameters of basis functions, one function per line. |

When calculations begin from scratch, only the header and command list need to be supplied.

Note that in any calculation that involves generation of a basis (either from scratch or by extending an existing basis) the input file also serves as an output file (hence the name - `inout.txt`). This file is updated with any incremental changes of the basis. The frequency with which the file gets overwritten can be controlled by the user.

In many practical situations one must save intermediate basis sets. This can be easily achieved with a command `SAVE_FILE`.

## Input file sections

### Header

### Command list

### History

### Basis functions

## Creating input files

When a large basis needs to be generated it requires multiple actions (e.g. basis enlargements `BASIS_ENL` and cyclic optimizations `OPT_CYCLE`) that need to be repeated routinely. It may be tedious to manually write the sequence of corresponding commands and their multiple arguments. For this reason it is convenient to create a simple utility that creates a desired list of commands, or, even better, prepares a ready-to-use input file `inout.txt` to initiate calculations from scratch.

Directory `ecgpack/utilities/input_file_manipulation` contains a Python script `ecg_input_file_generator.py` that can be invoked to generate an input file for a desired system with a long list of repeated actions. One can modify this script as needed.

## Restarting calculations that generate a basis
