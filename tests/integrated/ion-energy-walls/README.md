# Radial ion energy walls

The ion-energy component pairs explicitly with conservative `polarisation_drift`.
Auxiliary AC boundary settings belong to the drift section.

Six nonlinear 8×8 budgets cover `fixed`, `free_slip` and `no_slip` with both
Neumann and Dirichlet AC conditions through t=400. At every output, minimum
ion density and local thermal energy (`Eion` minus both kinetic reservoirs)
must exceed `1e-3` of their initial means. The Reynolds current advects full vorticity with MC and uses a central
metric remainder. Wall particle, momentum and energy weights are uniform in Z,
using products of mean species mass partitions and mean wall-face quantities.
Signed budgets retain Reynolds and chi exchange. Viscous stress work remains
a diagnostic; the polarisation energy flux carries viscous drift work. All
particle, energy and vorticity budgets use open-wall solver-relative limits.
Minimum checks exercise Reynolds wall exchange, chi exchange, viscous current
and wall motion.

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
currents are included. The chi solve excludes the viscous rate. Polarisation transport consequently
carries the viscous drift work; it retains the physical current faces. Vorticity
balances include signed Reynolds wall exchange in every wall mode.

Perturbed histories save 1600 intervals of 0.25 through t=400, because the
retained Reynolds wall current introduces a measurable trapezoidal cadence
error in vorticity as well as energy budgets. Flush/verification runs every
eight outputs, retaining every saved state with a two-unit disk-flush cadence.
Per-case and CTest time ceilings are 4800 and 14400 s. This budget matrix does
not establish wall accuracy. The enabled Python
`test_reynolds_wall_curvature_converges` regression checks 32/64/128 spatial
refinement and exposes a nonconvergent Reynolds boundary truncation error.
