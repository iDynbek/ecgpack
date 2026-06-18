# Description

This directory contains the sample wave-function files `wf_state0.txt` and `wf_state1.txt`, corresponding to the initial and final states of the helium atom in the infinite-nuclear-mass approximation, respectively.

Each wave function is expanded in a basis of 10 ECG basis functions.

No additional input file, such as `inout.txt`, is required to run the transition-matrix-element calculation.

Using double-precision arithmetic, the calculation takes approximately 1 second on a single CPU core when compiled with `gfortran` and executed on an AMD Ryzen 7 7800X3D processor.

## Final results

Atomic units are used throughout.

For this sample calculation, the transition energy is

$$
\Delta E = 0.779\,287\,387.
$$

The calculated length-form coordinate matrix element is

$$
\left\langle S^e \middle| z \middle| P^o \right\rangle
=
0.395\,665\,167.
$$

The calculated velocity-form momentum matrix element is purely imaginary and is given by

$$
\left\langle S^e \middle| P_z \middle| P^o \right\rangle
=
i\,0.314\,193\,144.
$$

The corresponding line strengths and oscillator strengths are evaluated separately during the post-processing step, for example using an Excel worksheet.

The length-form line strength is

$$
S^{\mathrm{L}} = 0.469\,652\,774,
$$

and the velocity-form line strength is

$$
S^{\mathrm{V}} = 0.296\,151\,995.
$$

The corresponding oscillator strengths are

$$
f^{\mathrm{L}} = 0.243\,996\,322
$$

and

$$
f^{\mathrm{V}} = 0.253\,352\,828.
$$

For comparison, a high-accuracy reference value obtained using Hylleraas-type wave functions is

$$
f^{\mathrm{L}}_{\mathrm{ref}} = 0.276\,1647.
$$

This reference value is reported in Table 12.11 of the *Handbook of Atomic, Molecular, and Optical Physics* (2023).

The disagreement between the present length- and velocity-form results, as well as their deviations from the reference value, is primarily attributable to the very small ECG basis used in this sample calculation, which contains only 10 ECG functions for each state.