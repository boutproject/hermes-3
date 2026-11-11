# Radial ion energy walls

Six nonlinear 8×8 budgets cover `fixed`, `free_slip` and `no_slip` with both
Neumann and Dirichlet AC conditions through t=400. At every output, minimum
ion density and local thermal energy (`Eion` minus both kinetic reservoirs)
must exceed `1e-3` of their initial means. The vorticity Reynolds term is
MC-upwinded and Lambda is central. The unpinned corner curl carries species,
momentum and energy through physical walls with zero net charge. The signed
budgets include this exchange, with Neumann face weights, and chi and viscous
wall exchanges. All particle and energy budgets use the existing open-wall
solver-relative limits. Minimum checks exercise the solenoidal wall exchange,
open chi exchange, viscous current and wall motion.

Dirichlet runs save diagnostics, including the common polarisation faces and
the no-slip chi DC midpoint check. A stripped run retaining standard `ion_chi`
must have complete wall coverage; removing chi then makes open AC coverage
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

Release validation completes all 13 executions. Both free-slip cases reach
`t=400`: minimum thermal fractions are 0.81304/0.81303 and minimum density
fractions are 0.93161/0.94632 (Neumann/Dirichlet). All budgets, reconstruction
comparisons and restart checks pass. Each solenoidal-energy minimum in
`criteria.json` is one tenth of that case's observed cumulative absolute
exchange fraction; the runner passes with these gates enabled.
