# Radial ion energy walls

Six nonlinear budgets cover `fixed`, `free_slip` and `no_slip` with both
Neumann and Dirichlet AC conditions. Pressure and vorticity perturbations
vanish at the X faces; curvature drives the finite-amplitude flow. Each 8×8 budget
runs to time 400 with 200 outputs, spanning several initial shear times. Minimum checks enforce perpendicular-energy
change (calibrated to `1e-5`, with observed changes above `1e-4`), absolute open chi particle/energy flux, free-slip viscous vorticity
flux and no-slip wall motion. Absolute face fluxes are summed before time
integration, so cancellation cannot hide an inactive test.

All Dirichlet runs save diagnostics. The no-slip chi DC midpoint is compared
with saved `phi_wall_rate` at each wall. Removing diagnostics from one retained
run must produce incomplete wall coverage. Separate analytic runs check rigid
rotation at 16/32 and uniform drift at 16; perturbed fixed/no-slip checkpoint
and restart runs compare with uninterrupted histories (13 executions total).

```sh
python tests/integrated/ion-energy-walls/runtest --exe build/hermes-3
```

Recovery, closed-current vorticity, fixed potentials and diagnostic comparisons
use `1e-10`. Free-slip vorticity uses `open_vorticity`: its cumulative viscous
wall flux has an output-sampling error, observed as `3.247e-8` / `2.593e-8` at
`rtol=1e-7` with the original operator. The perturbed cases use `rtol=1e-7`
and `atol=1e-10`; `rtol=1e-8` exceeded 600 s in a free-slip case. Analytic
profiles retain tighter tolerances. Signed budgets include the retained Reynolds
corner curl and its species/enthalpy face weights. Sampled energy exchanges use
`open_energy`; halving cadence raises the new fixed-Neumann residual from
`2.705e-9` to `1.042e-8`. Upper limits keep a `1e-10` floor.

The new free-slip evolution reaches density floors and nonphysical local energy
by time 240. Increasing auxiliary chi corrections from 2 to 8 and 20 did not
remove the failure; the configured correction count remains unchanged. Negative
integrated Reynolds heating precedes the failure. This is an unresolved acceptance
failure, not a budget limit to calibrate away. All physics checks remain enabled.
`criteria.json` records valid observations; isolated `result.json` files and the
PR validation description record timings and failures. No parallel current or
sheaths are present.
