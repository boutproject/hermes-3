1D flow with recycling at the sheath
====================================

Evolves two species, electrons and ions, with separate pressures
(temperatures) but the same flow (zero current). The 1D system
has a no-flow boundary on the lower Y boundary, and a sheath on
the upper Y boundary.

The simulation is driven by a uniform source of heat and particles,
as sources in the `Nd+`, `Pe` and `Pd+` equations.

A simple model for the sheath recycling is included:

 - Ions recycle as atoms.
 - Atoms ionise to create ions.
 - Atoms diffuse across the magnetic field, which appears in 1D
   equations as an enhanced parallel diffusion.
 - A non-uniform grid is used, which packs more points near the target


Solver
------

If using the CVODE (Sundials) solver, preconditioning is strongly recommended.

Without preconditioning:

    ./hermes-3 -d examples/1D-recycling/ solver:use_precon=false

    2.000e+03      60874       1.23e+02    96.7    0.0    0.3    0.0    2.9
    4.000e+03      30359       6.29e+01    96.7    0.0    0.3    0.0    2.9
    6.000e+03      30426       6.52e+01    96.8    0.0    0.3    0.0    2.9

With preconditioning:

    ./hermes-3 -d examples/1D-recycling/ solver:use_precon=true

    2.000e+03       5840       1.18e+01    87.9    0.0    0.2    0.3   11.6
    4.000e+03       2670       5.41e+00    87.5    0.0    0.2    0.5   11.8
    6.000e+03       2786       5.75e+00    87.6    0.0    0.2    0.4   11.7
