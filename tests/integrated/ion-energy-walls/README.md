# Radial ion energy walls

Six nonlinear 8×8 budgets cover `fixed`, `free_slip` and `no_slip` with both
Neumann and Dirichlet AC conditions through t=400. At every output, minimum
ion density and local thermal energy (`Eion` minus both kinetic reservoirs)
must exceed `1e-3` of their initial means. Other minimum checks enforce flow
activity, absolute open chi exchange, viscous wall current and wall motion.
Pinned Reynolds corners make every weighted normal Reynolds face exactly zero;
the signed budgets include only the remaining chi and viscous wall exchanges.

Dirichlet runs save diagnostics, including the no-slip chi DC midpoint check.
An in-memory stripped run must report incomplete coverage. Rigid rotation at
16/32, uniform drift at 16, and perturbed fixed/no-slip checkpoint/restart
comparisons bring the full runner to 13 executions.

```sh
python tests/integrated/ion-energy-walls/runtest --exe build/hermes-3
python tests/integrated/ion-energy-walls/runtest --exe build/hermes-3 --reynolds-reconstruction central
```

Recovery, diagnostic comparisons, closed-wall balances and fixed potentials
use `1e-10`. Free-slip vorticity uses the unchanged solver-relative
`open_vorticity` criterion for sampled viscous wall current. Dirichlet budgets
retain the existing open-energy/particle limits. No physical or solver
parameters are adjusted for the reconstruction comparison.

With pinning, XPPM free slip remains unstable: Neumann evolution stalls after
t=152 and Dirichlet potential recovery fails near t=218.17. Central free slip
completes both cases to t=400; its density/thermal minima are above the new
checks. Its interior Reynolds heating remains weakly negative, but much smaller
than XPPM. The default remains `xppm`; energy-consistent dissipation is a later
model decision. Central Neumann particles balance to roundoff, while sampled
viscous wall power leaves an energy residual of `2.40e-8`, failing the requested
`1e-10` floor. All acceptance checks remain enabled. `criteria.json` records
observations; the PR validation description gives timings, regional heating
and the pre-change wall-exchange diagnosis. Figures are outside the repository.
