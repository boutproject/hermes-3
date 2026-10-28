# Driven multispecies energy fluxes

This source-free X-Z slab evolves electron pressure with `evolve_pressure` and
total ion energy with `evolve_vorticity_ion_energy`. The two ions have
`(A,Z)=(2,1)` and `(3,2)`, unequal parallel momenta, spatially varying composition,
and a common temperature. `average_atomic_mass=1.75` differs from the actual
local mass density used in the current weights. Multiple independent X-Z modes
in density/pressure and an initial vorticity perturbation produce nonzero
Arakawa, Reynolds, electric and polarisation transport.

[data/BOUT.inp](data/BOUT.inp) uses 16 by 16 physical X-Z cells, four identical
Y planes and two guards in X/Y. X, Y and Z are periodic. The physical duration
is 120 normalised time units with output cadences 12, 6 and 3. Normalisations
are `Nnorm=2e18 m^-3`, `Tnorm=5 eV`, `Bnorm=0.35 T`.

## Power balance

The plasma inventory is saved ion total energy plus electron thermal energy.
External power and physical boundary power vanish. Electron/ion curvature work
redistributes energy internally. For the ion inventory, electron electrical
exchange `-phi*Div(Je)` is retained explicitly, together with the independently
reconstructed ExB, Reynolds, electric, parallel and polarisation face fluxes.
The ion residual is `delta Eion - integral(Eion_src-phi*Div(Je)) + integral Pout`.

The Python reconstruction uses MC transport and face-interpolated species
weights, the spectral `Delp2(phi)` in the Reynolds current, the Arakawa bracket,
and an independent sparse finite-volume solve for `ddt(Phi)`. It does not read
`ddt(Phi)` or other optional diagnostics to perform that solve. The periodic
potential gauge is removed only when comparing diagnostics. Both field and
volume-integrated errors are reported. All intended transport terms must exceed
fixed minimum amplitudes, and the three-cadence sampling differences must
shrink by at least a factor of two.

## Tracked discretisation limitations

Hermes XPPM suppresses the radial periodic seam with `bndry_flux=false` because
it tests `firstX/lastX`. Reconstruction matches that behaviour; the case checks
energy budgets without claiming correct transport across that seam. Reynolds
transport uses the same X-Z face operator as potential recovery and polarisation,
including Z fluxes when `MYG=0`. Flux reconstruction currently requires uniform
orthogonal slab metrics. The separate Kperp/forward operator tests include
nonuniform metrics and `g13`; full nonorthogonal transport is a later stage.

## Running and reports

From the repository root, use the workspace interpreter:

```sh
.venv/bin/python tests/integrated/2D-ion-energy-fluxes/runtest --exe build/hermes-3
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
