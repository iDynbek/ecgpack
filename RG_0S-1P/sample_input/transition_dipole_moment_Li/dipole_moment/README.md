# Description

This directory contains the sample wave-function files `wf_state0.txt` and `wf_state1.txt`, corresponding to the initial $2~^2S^e$ and final $2~^2P^o$ states of the lithium atom in the infinite-nuclear-mass approximation, respectively.

Each wave function is expanded in a basis of 100 ECG basis functions.

No additional input file, such as `inout.txt`, is required to run the transition-matrix-element calculation.

Using double-precision arithmetic, the calculation takes approximately 1 second on a single CPU core when compiled with `gfortran` and executed on an AMD Ryzen 7 7800X3D processor.

## Final results

Atomic units are used throughout.

For this sample calculation, the transition energy is

$$
\Delta E = 0.068\,044\,261.
$$

The calculated length-form coordinate matrix element is

$$
\left\langle S^e \middle| z \middle| P^o \right\rangle
=
2.342\,866\,004.
$$

The calculated velocity-form momentum matrix element is purely imaginary and is given by

$$
\left\langle S^e \middle| P_z \middle| P^o \right\rangle
=
0.159\,590\,277\,i.
$$

The corresponding line strengths and oscillator strengths are evaluated separately during the post-processing step, for example using an Excel worksheet.

The length-form line strength is

$$
S^{\mathrm{L}} = 32.934\,1266,
$$

and the velocity-form line strength is

$$
S^{\mathrm{V}} = 0.152\,814\,341.
$$

The corresponding oscillator strengths are

$$
f^{\mathrm{L}} = 0.746\,992\,775,
$$

and

$$
f^{\mathrm{V}} = 0.748\,602\,639.
$$

For comparison, high-accuracy nonrelativistic reference values for the infinite-nuclear-mass lithium atom, obtained using 11,000 and 12,000 ECG functions for the $2~^2S^e$ and $2~^2P^o$ states, respectively, are

$$
f_{\mathrm{ref}}^{\mathrm{L}} = 0.746\,956\,8098,
$$

and

$$
f_{\mathrm{ref}}^{\mathrm{V}} = 0.746\,956\,79.
$$

These reference values are reported in:

S. Nasiri, J. Liu, S. Bubin, M. Stanke, A. Kędziorski, and L. Adamowicz, *Atomic Data and Nuclear Data Tables* **149**, 101559 (2023), https://doi.org/10.1016/j.adt.2022.101559.

The disagreement between the present length- and velocity-form oscillator strengths, as well as their deviations from the reference values, is primarily attributable to incomplete basis-set convergence because only 100 ECG functions are used for each state in this sample calculation.
