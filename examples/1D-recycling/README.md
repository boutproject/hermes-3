1D flow with recycling at the sheath
====================================

Evolves two species, electrons and ions, with separate pressures
(temperatures) but the same flow (zero current). The 1D system
has a no-flow boundary on the lower Y boundary, and a sheath on
the upper Y boundary.

The simulation is driven by a uniform source of heat and particles,
as sources in the `Ni`, `Pe` and `Pi` equations.

A simple model for the sheath recycling is included:

 - Ions recycle as atoms.
 - Atoms ionise to create ions.
 - Atoms diffuse across the magnetic field, which appears in 1D
   equations as an enhanced parallel diffusion.
 - A non-uniform grid is used, which packs more points near the target

Some important limitations

 - No charge exchange between ions and atoms
 - No recombination
 - No molecules etc.

