# Open parallel-boundary energy transport

An all-open rectangular Y slab prescribes ion velocity `1000 m/s`, evolves its
density and total ion energy, and evolves electron pressure with quasineutral
electron density. The initial ion energy rises linearly along Y. Electron
parallel velocity is zero. `ixseps1=ixseps2=0` and the unbranched Y connections
identify a simple open rectangle, which the loader distinguishes from tokamak
geometry. Saved guards are retained; no sheath boundary values are fabricated.

[data/BOUT.inp](data/BOUT.inp) uses eight physical X cells, 16 Y cells, singleton
Z, periodic X, identity coordinates and homogeneous Neumann pressure/density
boundaries. The duration is 20 normalised time units, at cadences 2, 1 and 0.5.
The velocity is imposed to isolate transport; this is not a parallel force
balance, sheath or conduction test.

## Power balance

External sources vanish. Ion enthalpy plus kinetic energy uses the MC parallel
flux in the interior and the actual midpoint boundary rule at the two open
ends. Outward power is the high-Y integrated face power minus the low-Y power.
The electrical flux `phi*Ji` is also reconstructed and included. All supported
Reynolds, ExB and polarisation boundary powers are measured from their faces,
even when they vanish in this controlled input.

The plasma inventory change plus the integrated net outward face power must
vanish. The net outward energy must be positive and exceed the balance residual
by a factor of 100. Discrete Gauss checks compare face powers with independent
volume divergences, and radial subvolumes check shared-transfer cancellation.
Future tests can add open momentum, sheath heat losses and conduction after
those discrete boundary contracts are supported.

## Running and reports

From the repository root, use the workspace interpreter:

```sh
.venv/bin/python tests/integrated/2D-ion-energy-boundary/runtest --exe build/hermes-3
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
