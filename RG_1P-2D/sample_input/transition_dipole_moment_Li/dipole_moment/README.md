# Description

This directory contains the sample wave-function files `wf_state1.txt` and `wf_state2.txt`, corresponding to the initial $2~^2P^o$ and final $3~^2D^e$ states of the lithium atom in the infinite-nuclear-mass approximation, respectively.

Each wave function is expanded in a basis of 100 ECG basis functions.

No additional input file, such as `inout.txt`, is required to run the transition-matrix-element calculation.

Using double-precision arithmetic, the calculation takes approximately 1 second on a single CPU core when compiled with `gfortran` and executed on an AMD Ryzen 7 7800X3D processor.

## Final results

Atomic units are used throughout.

For this sample calculation, the transition energy is

$$
\Delta E = 0.074\,525\,543.
$$

The calculated length-form coordinate matrix element is

$$
\left\langle S^e \middle| z \middle| P^o \right\rangle
=
-2.219\,088\,686
.
$$

The calculated velocity-form momentum matrix element is purely imaginary and is given by

$$
\left\langle S^e \middle| P_z \middle| P^o \right\rangle
=-0.169\,427\,277\;i.
$$

The corresponding line strengths and oscillator strengths are evaluated separately during the post-processing step, for example using an Excel worksheet.

The length-form line strength is

$$
S^{\mathrm{L}} = 73.865\,318\,98,
$$

and the velocity-form line strength is

$$
S^{\mathrm{V}} = 0.430\,584\,032.
$$

The corresponding oscillator strengths are

$$
f^{\mathrm{L}} = 0.612\,507\,889,
$$

and

$$
f^{\mathrm{V}} =0.641\,064\,595.
$$

For comparison, high-accuracy nonrelativistic reference values for the infinite-nuclear-mass lithium atom, obtained using more than 32,000 Hylleraas-type basis functions for the $2~^2P^o$ and $3~^2D^e$ states, respectively, are

$$
f_{\mathrm{ref}}^{\mathrm{L}} = 0.638\,568\,129
,
$$

The reference value is reported in:

Deng Sun , Liming Wang , Zong-Chao Yan, *Atomic Data and Nuclear Data Tables* **149**, 101559 (2023), https://doi.org/10.1016/j.adt.2026.101793.

The disagreement between the present length- and velocity-form oscillator strengths, as well as their deviations from the reference values, is primarily attributable to incomplete basis-set convergence because only 100 ECG functions are used for each state in this sample calculation.
