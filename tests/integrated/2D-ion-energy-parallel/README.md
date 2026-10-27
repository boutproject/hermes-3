# Electron parallel energy and a tracked heating defect

A periodic rectangular slab evolves electron density, parallel momentum and
pressure. A single quasineutral ion has stationary total ion energy, uniform
ion pressure and the common recovered temperature. The electron initial
momentum is `0.0025*sin(y)`; density and electron pressure start at 1.
[data/BOUT.inp](data/BOUT.inp) uses 16 Y cells, eight physical X cells, singleton
Z, identity parallel coordinates, MC transport and numerical viscous heating.
The duration is 2 normalised time units, at cadences 0.2, 0.1 and 0.05.

## Energy accounting

The plasma inventory includes evolved electron parallel kinetic energy as well
as electron thermal and ion total energy. The runner independently reconstructs
parallel split fluxes, density and momentum RHS, `v*Grad_par(Pe)`, and numerical
viscous heating. The kinetic RHS uses the actual evolved electron density
derivative: `v*ddt(NVe)-A*v^2*ddt(Ne)/2`. Pressure work and kinetic energy must
be measurably nonzero. Electron thermal balance, ion balance and all available
RHS diagnostics are checked independently of the saved diagnostics.

This input isolates parallel pressure work and numerical transport; it does not
exercise parallel electrical work or ion sound response. Electrical X-Z work
is exercised in the separate driven-flux test. Open electron momentum,
sheath and conduction budgets are later additions.

## Tracked BOUT numerical-heating discrepancy

The present `FV::Div_par_fvv_heating` uses `n_edge*(a +/- v_edge/2)`, whereas the
split particle/momentum flux requires `n_edge*(a +/- v_edge)/2` in the matching
kinetic balance. Its excess is `n_edge*a/2` times the velocity-difference factor.
A synthetic operator test establishes this identity, and the integrated test
independently predicts the resulting excess plasma power.

The report preserves the raw external-source/face plasma balance and explicitly
marks `external_plasma_balance_conserves=false`. It separately reports the
excess-power history and the residual after accounting for that predicted
numerical defect. The latter must be below `1e-8` of reference energy, and the
raw drift must exceed it by a factor of 20. Thus a passing test means the defect
is quantified and diagnostics agree; it does **not** claim source-free plasma
conservation. At the finest checked cadence, the raw drift is approximately `1.43e-6` of
reference energy, while the independently defect-accounted residual is about
`2e-12`. No BOUT physics implementation is changed by this test. Removing
or fixing the defect requires updating this named regression explicitly.

## Running and reports

From the repository root, use the workspace interpreter:

```sh
.venv/bin/python tests/integrated/2D-ion-energy-parallel/runtest --exe build/hermes-3
```

The identically named CTest runs from the configured source/build test directory.
`--profile ci` (default) runs three output cadences at fixed duration and a
finest-cadence diagnostics-off pair. `--profile smoke` uses two cadences and
the off pair. `--output-root DIR` selects an existing results parent. Every
invocation creates a fresh directory and retains input hashes, execution status,
logs, per-run `energy-report.json` and a final `result.json`. MPI requires local
network transport even with one rank; this Mac checkout uses
`FI_SOCKETS_IFACE=en0`.

[criteria.json](criteria.json) fixes inputs, SI source rates, exercised signals
and acceptance limits. The shared implementation is
[ion_energy_flux_test.py](../ion_energy_flux_test.py); the independent Python
operators are in [energy_analysis](../../../tools/energy_analysis/README.md).
State and inventory histories must agree between diagnostic pairs. Removing
optional time-dependent fields, while preserving saved closure momenta and
prescribed velocities, must leave the reconstructed fluxes and inventories
unchanged. Sources are established from explicit input rates and a component
whitelist, never from missing diagnostics.

All records, including the initial record, check the ion partition and forward
operator (`1e-9` normalised absolute limits), flux/RHS diagnostic recovery
(`1e-10`), and discrete Gauss balance (`1e-10` volume-integrated normalised
limit). Species charge face fluxes and divergences must sum to the common
polarisation current. Two radial subvolumes independently check ion balances
and cancellation of their shared transfer. Cumulative and interval ion/plasma
balances use `1e-6` relative limits; regional balances use `2e-6`. Reports keep
sources, outward face powers, exchanges, sampling estimates and unavailable
coverage separate. These cases add no production diagnostic arrays.
