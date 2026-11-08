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

Recovery, vorticity, fixed wall potentials and diagnostic state comparisons
use `1e-10`. Open budgets, restart and analytic rates use limits resolved from
each run's `solver:rtol`, with a `1e-10` floor. The analytic profiles retain
tighter solver tolerances. `criteria.json` records observations; isolated
`result.json` files contain effective limits and times. No parallel current
or sheaths are present. Physics failures remain enabled.

Known failures: free-slip Neumann/Dirichlet vorticity residuals are
`3.244e-8` / `2.593e-8` against `1e-10`. The cumulative snapshot-flux checks
remain enabled; all magnitude, chi-hook, analytic and restart checks pass.
