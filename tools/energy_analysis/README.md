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

`zero_current_defect` reports the unpaired electron parallel force work
`Ve Grad_par(Pe) + Je Grad_par(phi)` in watts. A `zero_current` electron closure
sets velocity from ion current without electron momentum balance or an evolved
kinetic reservoir; zero total current does not imply zero work. The supported
reconstruction uses the default `v*Grad(p)` pressure form, C2 parallel gradients
and uniform Y metrics. It requires retained Y guards, or explicit
`metadata["periodic_y"] = True` for a guard-free periodic dataset. An exactly
Y-invariant pressure and potential give zero work without reconstruction.
`closed_balance` exposes `model_defect_power`; `source_balance` retains
`known_model_defect_power` separately from external heating. Nonzero or
unsupported closure work makes these balances incomplete and leaves a closed
source residual unavailable. The strict X-Z closed test uses Y-invariant
fields; its earlier Y-varying run is retained as evidence of this limitation.

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
- Polarisation particle, momentum and energy wall exchanges from saved `phi`,
  `Vort`, species state and the standard `ion_chi` output; stored wall midpoint
  histories and viscous wall fluxes.
- Electron curvature work for explicit constant Z curvature and FFT derivatives.

`ion_chi` is saved even with `diagnose=false`. `wall_budget` reconstructs
full-vorticity MC transport plus the central Reynolds remainder. It applies
Z-uniform wall weights: adjacent-cell mean mass partitions, separate means of
wall-face temperature and parallel velocity, and the square of mean velocity
for parallel kinetic energy. It reports signed exchanges by wall and physical
Y row, their totals and absolute face-flux sums. The net work is the wall
midpoint of `phi+psi` times wall current. Full-current budgets require `ion_chi`
with both AC conditions;
snapshots without it report incomplete coverage for that contribution. Sheaths,
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

`reconstruct_fluxes` rebuilds full-vorticity MC advection and the central
metric remainder, then combines Reynolds and chi faces into the polarisation
current. Particle, momentum and energy use arithmetic interior face averages
and the uniform wall products described above. The energy coefficient is
`phi+psi+sum_s[(A_s/Z_s)/rho*(1.5*P_s+Kparallel_s)]`.
The optional `polarisation_flux_x/z` diagnostics include area integration and
use the lower-cell face index; their comparisons include the lower X wall and
exclude undefined outer guards. Saved chi is used directly; an auxiliary solve
supports older snapshots that lack it.

`DivEExB` advects `Eion`. Post-processing provides
`DivEExB_perpendicular`, `DivEExB_thermal_parallel` and `reynolds_heating`.
The perpendicular part uses MC reconstruction of the discrete Kperp. Since MC
is nonlinear, the thermal-plus-parallel part is the total flux minus this part,
so the two add exactly, including their reconstruction errors. Heating is
`reynolds_work - DivEExB_perpendicular`; only `reynolds_work` is saved by the
solver. Full flux reconstruction retains the uniform orthogonal slab contract;
the face-based perpendicular operator itself supports varying J, spacing and
g13. The wall budget recovers full Reynolds/chi wall currents on rectangular
grids with varying metrics without requiring complete RHS coverage.

## Component version and radial currents

The output tag `evolve_vorticity_ion_energy_version` identifies the equation
and discretisation contract. Version 5 moves conservative transport to
`polarisation_drift` version 1, selected explicitly with
`implementation=conservative`. Its `polarisation_ac_boundary`,
`polarisation_corrections` and `laplacian` settings belong to the drift section.
Version 4 and untagged snapshots retain their original ion-energy ownership.
Version 6 pairs with conservative drift version 2, includes completed parallel
current in auxiliary/vorticity evolution and evaluates parallel electrical
work after sheath conditioning. Versions 4/5/6 are supported. The existing
slab budget coverage remains restricted; version support does not establish
tokamak/sheath power-ledger coverage. Unknown versions are rejected. Reynolds and chi wall faces are retained for
both AC conditions. Numerical diagnostic names and conversions are preserved.

The auxiliary solve excludes the viscous rate, so polarisation transport
carries ion viscous drift work. `DivEviscous` remains diagnostic stress work;
it is excluded from the energy RHS and outward energy budget to avoid double
counting. Electron heating remains an internal transfer. Saved chi takes
precedence over reconstruction. The default is ten defect corrections.

For versions 4/5, current continuity is `DivJpol+DivJdia=Rchi`, including wall
cells. Version 6 adds `DivJpar` and enabled `DivJextra`/`DivJcol` on the left.
Its existing X-Z analysis paths are tested; general nonzero-parallel-current
RHS/chi reconstruction remains unfinished and is not validated by tag support. The residual
controls the gauge dependence of local energy partition. Global budgets alone
do not establish radial boundary truncation accuracy. Chi's DC wall-rate
approximation still neglects the mean pressure-potential rate.

Use `PYTHONPATH=tools python -m energy_analysis.spectra RUN_DIR --output bands.json`
for fixed-mode wall/interior RMS, radial profiles and e-folding widths. This
streaming monitor supports singleton-Y, X-decomposed native version-4/5/6 dumps;
supply physical domain lengths if they differ from 50 mm. Bands 16–31 and
32–63 are fixed physical modes, independent of resolution. Wall amplitudes use
first physical cell centres. Widths include their half-cell offset from the
wall; crossings without wall dominance are flagged. A missing crossing has
no width. Native pytest regressions run with `HERMES_EXECUTABLE=build/hermes-3`;
add `--run-slow` for the drift-plane envelope test.
