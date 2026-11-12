# Radial ion energy walls

Six nonlinear 8×8 budgets cover `fixed`, `free_slip` and `no_slip` with both
Neumann and Dirichlet AC conditions through t=400. At every output, minimum
ion density and local thermal energy (`Eion` minus both kinetic reservoirs)
must exceed `1e-3` of their initial means. The vorticity Reynolds term is
MC-upwinded and Lambda is central. The unpinned corner curl carries species,
momentum and energy through physical walls with zero net charge. The signed
budgets include this exchange, with Neumann face weights, and chi and viscous
wall exchanges. All particle, energy and vorticity budgets use the existing open-wall
solver-relative limits. Minimum checks exercise the solenoidal wall exchange,
open chi exchange, viscous current and wall motion.

Dirichlet runs save diagnostics, including the common polarisation faces and
the no-slip chi DC midpoint check. A stripped run retaining standard `ion_chi`
must have complete wall coverage; removing chi then makes both AC conditions' coverage
incomplete. Rigid rotation at 16/32, uniform drift at 16, and perturbed
fixed/no-slip checkpoint/restart comparisons bring the runner to 13 executions.

```sh
python tests/integrated/ion-energy-walls/runtest --exe build/hermes-3
```

Recovery, diagnostic comparisons and fixed potentials use `1e-10`. Free-slip
vorticity uses the unchanged solver-relative `open_vorticity` limit for sampled
viscous wall current. Solver parameters, density/thermal minimum checks and
existing tolerances remain unchanged. Free-slip execution failures are reported
with saved progress; they are not addressed by tuning inputs or thresholds.
Validation timings and observed wall exchanges belong in the PR description.

The earlier validation, before the default-off `include_divjviscous` option,
completed all 13 executions. These recorded values use viscosity in chi and
do not validate the current default. Both free-slip cases reached
`t=400`: minimum thermal fractions are 0.81304/0.81303 and minimum density
fractions were 0.93161/0.94632 (Neumann/Dirichlet). All budgets, reconstruction
comparisons and restart checks passed. Each solenoidal-energy minimum in
`criteria.json` is one tenth of that case's observed cumulative absolute
exchange fraction; the runner passes with these gates enabled.

The complete Reynolds and solved chi wall currents are retained for both AC
conditions. The AC choice controls chi potential guards; DC and cross-metric
currents are included. Q retains odd wall guards and Lambda linear extrapolation. Vorticity balances include signed Reynolds wall
exchange, so all modes use the existing open-vorticity sampling limit. The
curl-only zero-charge and weighted exchange checks remain separate.

Perturbed histories now save 1600 intervals of 0.25 through t=400, because the
retained Reynolds wall current introduces a measurable trapezoidal cadence
error in vorticity as well as energy budgets. Solver tolerances and budget
limits are unchanged. Flush/verification runs every eight outputs, retaining
every saved state and the previous two-unit disk-flush cadence. Per-case and
CTest time ceilings are scaled by the eightfold increase in output count to
4800 and 14400 s. This budget matrix does not establish wall regularity;
the active manufactured boundary-order regression and matched 32/64/128
refinement expose a remaining viscous wall defect.

The full-current validation completed all 13 executions and passed all 204
checks. Maximum nonlinear energy, particle and vorticity budget residuals are
5.17e-9, 2.39e-9 and 6.17e-9 per volume, respectively. Both checkpoint/restart
comparisons pass. These budget results are separate from the failed wall-order
and viscous-refinement acceptance checks.
