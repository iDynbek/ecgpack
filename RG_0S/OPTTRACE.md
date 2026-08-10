# Optimizer trajectory tracing (for animation)

Records what the nonlinear optimizer actually *does* — every candidate parameter
set it tries, and the energy it got — so the optimization of a basis function can
be replayed as an animation.

Off by default and free when off. This is instrumentation: it lives on the
`viz-opttrace` branch and is deliberately not on the production branches.

## Use

```bash
export ECG_OPTTRACE=1                    # enable
export ECG_OPTTRACE_FILE=mytrace.csv     # optional, default opttrace.csv
mpirun -np 1 ./ecg
```

Only rank 0 writes, so the file is the same regardless of rank count (though the
*trajectory* itself can differ slightly between rank counts — floating-point
summation order changes, and the optimizer's accept/reject decisions amplify it).

Currently hooked into `OptCycleI` and `OptCycleG` — the routines that refine
already-existing basis functions. Growth (`BasisEnl*`) and `FullOpt1*` are not
traced yet; the same two-line pattern would extend it.

## Output

Two files. `opttrace.csv` has one row per objective evaluation per optimized
function:

| column | meaning |
|---|---|
| `step` | global counter, one per objective evaluation |
| `phase` | which routine produced the row (`OPTCYCLE_I` / `OPTCYCLE_G`) |
| `cbs` | basis size at the time |
| `slot` | index within the set being optimized together (1..nfo) |
| `func` | the function's own number in the basis |
| `kind` | 1 = optimizer asked for an energy, 2 = it asked for a gradient |
| `status` | 0 = evaluation succeeded; nonzero = failed, `energy` is meaningless |
| `energy` | variational energy at this point |
| `p1..pN` | the function's nonlinear parameters at this point |

`opttrace.csv.meta` carries what you need to interpret it: particle count,
masses, charges, and the map from `p1..pN` to the matrix `L`.

**Filter on `status == 0`.** Failed evaluations are recorded deliberately (a
failed inverse iteration is part of the story) but their energies are garbage.

## Turning parameters into something you can draw

Each parameter vector is the lower triangle of `L`, stored column by column:

```
indx = 0
for i in 1..n:
    for j in i..n:
        indx += 1
        L[j][i] = p[indx]
```

The Gaussian's exponent matrix is `A = L @ L.T` (symmetric, positive definite by
construction — that's *why* the code optimizes `L` rather than `A` directly), and
the basis function is

```
f(r) = exp( -r^T A r )
```

where `r` stacks the pseudoparticle coordinates: `r_i = x_(i+1) - x_1`, i.e. each
particle's position relative to particle 1 (the nucleus, for an atom).

### The easy case to animate: helium

He has 3 particles → 2 pseudoparticles → `A` is 2×2 → 3 parameters. The exponent
is

```
A11*|r1|^2  +  2*A21*(r1 . r2)  +  A22*|r2|^2
```

`r1` and `r2` are the two electron positions relative to the nucleus. That's still
6-dimensional, so pick a slice to draw. Two that look good:

- **Fix the opening angle** between the electrons (say 90°, so `r1 . r2 = 0`) and
  heat-map the density over `(|r1|, |r2|)`. Each animation frame is one optimizer
  step; you watch the blob stretch and rotate as `A` changes.
- **Plot the three parameters as a path** in 3-D, coloured by energy. This shows
  the optimizer's search directly — including the failed excursions.

The off-diagonal `A21` is the interesting one physically: it is the *correlation*
between the two electrons. A diagonal `A` would be an uncorrelated product of two
independent Gaussians; `A21 != 0` is exactly what makes these "explicitly
correlated" Gaussians worth the trouble.

`sample_input/basis_generation_He_1Se-01` is a good source case — small, fast, and
its exact energy (−2.9037) is known, so the descent has a meaningful target.
