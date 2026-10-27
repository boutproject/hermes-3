# Prescribed electron thermal heating

A uniform two-ion slab uses common-temperature total ion energy and electron
`evolve_pressure`. The ion masses/charges are `(2,1)` and `(3,2)`, with densities
0.5 and 0.25 and opposite parallel momenta. All fields are uniform, so the
periodic face powers and electron/ion curvature exchange vanish despite the
configured curvature vector.

[data/BOUT.inp](data/BOUT.inp) prescribes `Pe:source=1e5 Pa/s`. Pressure is
`2/3` of thermal energy density, so the deposited power density is exactly
`1.5e5 W/m^3`. No ion, particle or momentum source is present. The source is
explicitly validated against the input; no diagnostic is needed to recover it.
The run lasts 120 normalised time units, with cadences 12, 6 and 3.

## Power balance

The plasma inventory must increase by `1.5 * Pe:source * simulated_volume * t`,
while saved ion energy remains constant. The runner also checks the analytic
electron thermal increase and requires deposited energy to exceed the measured
balance residual by a factor of 100. `Eion_src`, `ddt(Pe)` and the energy-flux
fields cross-check the independent calculation when diagnostics are enabled.
The complete balance still works in the off and stripped-output cases.

A constant source is integrated exactly by the trapezoidal rule. Cadence
reports therefore quantify solver/roundoff error without requiring a fictitious
fourfold reduction of a zero sampling error.

## Running and reports

From the repository root, use the workspace interpreter:

```sh
.venv/bin/python tests/integrated/2D-ion-energy-heating/runtest --exe build/hermes-3
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
