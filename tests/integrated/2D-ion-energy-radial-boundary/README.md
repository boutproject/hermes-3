# Radial wall models

`runtest` checks every `radial_boundary` × `polarisation_ac_boundary` combination
in a seeded, orthogonal x–z slab with no sources, curvature or parallel dynamics.
The free-slip cases include viscosity. Particle content is checked against
`particle_pol_wall_loss_h+`, total energy against the open polarisation and
viscous face fluxes, and vorticity against viscous diffusion. It also checks
potential and vorticity guards and a native restart round trip in every mode.

The no-slip rigid-rotation cases use 16, 32 and 64 physical radial cells.
A linear vorticity profile excites a cubic potential plus a relaxing linear
component. Cold, uniform ions suppress the change in the mean pressure
potential from enthalpy transport, isolating the prescribed E×B wall law. The fitted time must match `L_phys/(2*kappa)` within 1%, with a
resolution spread below 0.5%. These cases retain the time solver and native
wall storage rather than integrating a separate wall ODE.

Run with `.venv/bin/python tests/integrated/2D-ion-energy-radial-boundary/runtest
--exe build/hermes-3 --output-root /private/tmp`. Each run writes its inputs,
outputs and logs, and the study writes `result.json` after all checks pass.

The existing quadratic-wall regression checks the second-order wall stencil,
and the tokamak relaxation test checks automatic closed-core averaging.
The independent slab flux reconstruction covers the selected AC conditions.
Parallel-current terms and the alternative Reynolds split remain outside Plan C.
