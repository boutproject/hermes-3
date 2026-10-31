# Two-ion MMS with kinematic viscosity

This focused variant of [the slab MMS](../3D-ion-energy-mms/README.md) uses
identity metrics and normalised kinematic viscosity `nu=0.02`. With the
chosen `Tnorm=Bnorm=1`, its SI input value is also `0.02 m^2/s`.
X remains nonperiodic; Y and Z remain periodic, without sheaths.

The exact potential is zero. The pressure potential varies as
`(1-cos(2*pi*X/10))^3`, so vorticity and its normal derivative vanish at
both X faces. The analytic current-continuity potential is
`chi=nu*Vort/2`: it satisfies the auxiliary midpoint Dirichlet DC / radial
Neumann AC boundary rule. The polarisation current has divergence
`-nu*lap(Vort)`, which is retained in all weighted density, momentum and energy
fluxes.

Sources include the viscous total-energy flux and the electron share of
viscous heating. Ion heating enters through the kinetic/thermal partition;
it is not added again to total ion energy. Numerical parallel viscosity and
its heating remain active without adding discrete cancellation sources.

Generate all inputs with `.venv/bin/python
tests/integrated/3D-ion-energy-mms/mms.py --all --check`, then run
`.venv/bin/python tests/integrated/3D-ion-energy-mms-viscosity/runtest` from
the repository root. The shared runner and analysis provide the same
resolution, tolerance, MPI, boundary-error and recovered-field checks as
the baseline.
