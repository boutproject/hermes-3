# Energy analysis

This is the shared analysis library for total ion energy tests and drift-plane
examples. It reads native BOUT++ rank dumps with `boutdata` and `xarray`, retains
physical guards, and uses resolved `BOUT.settings` (falling back to `BOUT.inp`).
An already annotated xhermes dataset can also be passed to `open_run`; xhermes
is not required when loading a directory.

```sh
python -m pip install -r tools/energy_analysis/requirements.txt
PYTHONPATH=tools python -m energy_analysis RUN_DIR --json report.json
```

```python
from energy_analysis import (
    open_run,
    energy_inventory,
    closed_balance,
    external_power,
    source_balance,
    wall_budget,
    compare_diagnostics,
    summary,
)

run = open_run("RUN_DIR")
inventory = energy_inventory(run)  # joules; physical time in seconds
balance = source_balance(run)  # cumulative energy and particle residuals
report = summary(run)  # JSON-serialisable digest
```

`energy_inventory` reconstructs ion thermal, parallel and perpendicular kinetic
energy, electron thermal energy, and electron parallel kinetic energy when its
momentum is evolved. It also exposes the saved total ion reservoir.
`closed_balance` reports energy, species particle and charge drifts from the
first record. `external_power` reconstructs pressure-source heating and the
moving-ion kinetic conversion of density sources, plus particle rates.
`source_balance` subtracts time-integrated external power and adds outward wall
power. `wall_budget` reconstructs viscous and polarisation wall fluxes and
measures stored midpoint-potential rates between output records. `summary`
includes inventories, particle content, recovery consistency, internal
curvature exchange, diagnostic comparisons and balance coverage.

All dataset results have `coverage` (`complete` or `incomplete`) and
`missing_terms`. Incomplete source balances retain known powers and inventory
changes but omit a residual. Unknown terms are never assigned a zero power.
`cadence_comparison` estimates trapezoidal sampling error on nested regular
output histories; its estimate assumes smooth powers. `AuditError` identifies
missing state or an unsupported reconstruction contract.

## Analysis without diagnostics

Complete from ordinary state, geometry and resolved options:

- Energy inventories and recovery consistency in supported slab geometry.
- Closed energy/particle/charge drift, including periodic domains.
- Static physical pressure and density source powers, with their SI conversion.
- Closed polarisation wall faces, stored wall midpoint histories, and viscous
  wall fluxes from saved `phi`, `Vort` and viscosity.
- Electron curvature work for explicit constant Z curvature and FFT derivatives.

Open AC polarisation wall faces need chi, which is normally a diagnostic.
When it is absent, that wall term reports incomplete coverage. Sheaths,
unrecognised source components, unsaved mesh sources, time-dependent source
prefactors, core-only source masks, and unimplemented geometry also report
missing terms. Collisional heat exchange cancels in the plasma total.
Electron friction under a zero-current closure requires mechanical work that
is not carried by an evolved electron reservoir and is not assumed to cancel.

Optional diagnostics are compared with reconstructed inventories. Existing
`fluxes`, `operators`, `transport`, `pressure` and `convergence` APIs remain
available for supported snapshot audits and the publication example.
`testing.Study` supplies isolated run directories, hashes, logs, execution status,
criteria checks and `result.json` to the thin integrated runners.

Run Python checks with `python -m pytest tests/python`, or use CTest's
`energy-analysis-python`. Plotting is optional; the
[publication example](../../examples/tokamak-2D-driftplane/2D-energy/README.md)
uses `matplotlib`.
