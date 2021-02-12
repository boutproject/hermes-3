1D flow with sheath boundary and heat conduction
================================================

Evolves two species, electrons and ions, with separate pressures
(temperatures) but the same flow (zero current). The 1D system
has a no-flow boundary on the lower Y boundary, and a sheath on
the upper Y boundary.

The simulation is driven by a uniform source of heat and particles,
as sources in the `Ni`, `Pe` and `Pi` equations.

This simulation does include heat conduction, and so also includes
electron-ion collisions which transfer energy between them, but
is otherwise the same as the `1D-sheath` example.
