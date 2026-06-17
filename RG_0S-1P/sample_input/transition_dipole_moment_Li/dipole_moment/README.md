# Description

This directory contains sample wave function files `wf_state0.txt` and `wf_state1.txt`, corresponding to the initial and final states of the Li atom (assuming infinite nuclear mass), respectively.
Each file contains a basis of 100 ECG functions.
To run the transition dipole moment calculation, no input file is required.
The calculation should take approximately 1–2 seconds on a single CPU core when using double precision (gfortran compiler on an AMD Ryzen 7 7800X3D).


## Final results

After running the `RG_0S-1P` code, the calculated transition quantities can be reported as follows:

$$
\Delta E = \text{value}
$$

$$
\left\langle S^e \middle| z \middle| P^o \right\rangle
=
\text{value}
$$

$$
\left\langle S^e \middle| P_z \middle| P^o \right\rangle
=
\text{value}
$$

The line strength and oscillator strengths are calculated separately in the post-processing step, for example in an Excel worksheet:

$$
S^{\mathrm{L}} = \text{value}
$$

$$
S^{\mathrm{V}} = \text{value}
$$

$$
f^{\mathrm{L}} = \text{value}
$$

$$
f^{\mathrm{V}} = \text{value}
$$