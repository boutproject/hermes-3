# Particle, momentum and thermal source conversion

A spatially uniform common-temperature mixture has two moving ions with
`(A,Z)=(2,1)` and `(3,2)`. Electrons use quasineutral density, zero-current
velocity and `evolve_pressure`. The test-only `test_energy_sources` component
adds independent constant density, momentum and thermal energy rates. It is
registered only when `HERMES_TESTS` is enabled and retains these values as scalar
options in `BOUT.settings`; it adds no time-dependent volume output.

[data/BOUT.inp](data/BOUT.inp) prescribes:

| Ion | Particle rate [m^-3 s^-1] | Momentum rate [N m^-3] | Thermal rate [W m^-3] |
| --- | --- | --- | --- |
| d+ | 6.7e22 | 1 | 1e4 |
| impurity | -3.35e22 | -0.5 | 5e3 |

The charge-weighted particle rates cancel, keeping electron density fixed
while the imposed momentum rates change velocities and the zero-current
electron closure. Actual ion mass density changes. The run
lasts 120 normalised time units, at output cadences 6, 3 and 1.5. Tight time-solver
settings resolve the independent cadence study.

## Power balance

An external `energy_source` is thermal. Its contribution to total ion energy is

$$Q_i=\sum_s\left[Q_{\mathrm{thermal},s}
  +v_s S_{\mathrm{momentum},s}
  -\frac12 A_s v_s^2 S_{\mathrm{density},s}\right].$$

Python reconstructs the rates from the manifest and saved normalisations,
and cross-checks `Eion_src` independently. Internal polarisation density and
momentum sources belong to the flux budget and are excluded from external
source conversion. Uniformity makes these internal sources zero here; the
separate driven-flux case exercises their nonzero current partition.

The analytic solution is `N=N0+Sn*t`, `NV=NV0+Sm*t`. Total ion energy changes by
known thermal deposition plus the change in summed parallel kinetic energy.
The runner checks this analytic energy change and each density history as well
as cumulative, interval and regional power balances. Output refinement checks
sampling differences separately from a persistent time-solver residual.

## Running and reports

From the repository root, use the workspace interpreter:

```sh
.venv/bin/python tests/integrated/2D-ion-energy-sources/runtest --exe build/hermes-3
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
