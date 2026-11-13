# Radial ion energy walls

Six nonlinear 8×8 budgets cover `fixed`, `free_slip` and `no_slip` with both
Neumann and Dirichlet AC conditions through t=400. At every output, minimum
ion density and local thermal energy (`Eion` minus both kinetic reservoirs)
must exceed `1e-3` of their initial means. The vorticity Reynolds term is
MC-upwinded and Lambda uses MC reconstruction with the velocity from
`Phi=phi+psi`. The unpinned corner curl carries species,
momentum and energy through physical walls with zero net charge. The signed
budgets include this exchange, with Neumann face weights, and chi and viscous
wall exchanges. All particle, energy and vorticity budgets use the open-wall
solver-relative limits. Minimum checks exercise the solenoidal wall exchange,
open chi exchange, viscous current and wall motion.

Dirichlet runs save diagnostics, including the common polarisation faces and
the no-slip chi DC midpoint check. A stripped run retaining standard `ion_chi`
must have complete wall coverage; removing chi makes coverage incomplete for
both AC conditions. Rigid rotation at 16/32, uniform drift at 16, and perturbed
fixed/no-slip checkpoint/restart comparisons bring the runner to 13 executions.

```sh
python tests/integrated/ion-energy-walls/runtest --exe build/hermes-3
```

Recovery, diagnostic comparisons and fixed potentials use `1e-10`. Free-slip
vorticity uses the solver-relative `open_vorticity` limit for sampled viscous
wall current. Free-slip execution failures are reported with saved progress.
Each invocation records wall exchanges and check results in `result.json`.

The complete Reynolds and solved chi wall currents are retained for both AC
conditions. The AC choice controls chi potential guards; DC and cross-metric
currents are included. Viscous stress is excluded from the chi solve and
polarisation transport.
Q retains odd wall guards and Lambda linear extrapolation. Vorticity balances
include signed Reynolds wall exchange, so all modes use the open-vorticity
sampling limit. The curl-only zero-charge and weighted exchange checks remain
separate.

Perturbed histories save 1600 intervals of 0.25 through t=400, because the
retained Reynolds wall current introduces a measurable trapezoidal cadence
error in vorticity as well as energy budgets. Flush/verification runs every
eight outputs, retaining every saved state with a two-unit disk-flush cadence.
Per-case and CTest time ceilings are 4800 and 14400 s. This budget matrix does
not establish wall accuracy. The enabled Python
`test_reynolds_wall_curvature_converges` regression checks 32/64/128 spatial
refinement and exposes a nonconvergent Reynolds boundary truncation error.
