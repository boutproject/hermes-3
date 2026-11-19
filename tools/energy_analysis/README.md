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
local wall weights: guard/cell midpoints of the mass partition, temperature and
parallel velocity; parallel kinetic energy squares
the local face velocity. Neumann pressure guards retain wall fluctuations while
the electrostatic midpoint is conducting and Z-uniform. The local work weight
is `Phi_w = phi_b + psi_w(z)`, so net work is the face sum of `Phi_w * Jpol_n`.
Signed exchanges are reported by wall and physical Y row, with totals and
absolute face-flux sums. Full-current budgets require `ion_chi` with both AC
conditions; snapshots without it report incomplete coverage for that contribution. Sheaths,
unrecognised source components, unsaved mesh sources, time-dependent source
prefactors, core-only source masks, and unimplemented geometry also report
missing terms. Collisional heat exchange cancels in the plasma total.
Electron friction under a zero-current closure requires mechanical work that
is not carried by an evolved electron reservoir and is not assumed to cancel.

Optional diagnostics are compared with reconstructed inventories. Existing
`fluxes`, `operators`, `transport`, `pressure` and `convergence` APIs remain
available for supported snapshot audits and the resolution-scan example.
`testing.Study` supplies isolated run directories, hashes, logs, execution status,
criteria checks and `result.json` to the thin integrated runners.

Run Python checks with `python -m pytest tests/python`, or use CTest's
`energy-analysis-python`. Plotting is optional; the
[resolution-scan example](../../examples/tokamak-2D-driftplane/2D-energy/README.md)
uses `matplotlib`.

`reconstruct_fluxes` rebuilds full-vorticity MC advection and the central
metric remainder, then combines Reynolds and chi faces into the polarisation
current. Particle, momentum and energy use arithmetic interior face averages
and the local wall weights described above. The energy coefficient is
`phi+psi+sum_s[(A_s/Z_s)/rho*(1.5*P_s+Kparallel_s)]`.
The optional `polarisation_flux_x/z` diagnostics include area integration and
use the lower-cell face index; their comparisons include the lower X wall and
exclude undefined outer guards. Saved chi is used directly; an auxiliary solve
is available for supported version-1 snapshots without chi within the slab
RHS reconstruction coverage.

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

The output tag `evolve_vorticity_ion_energy_version = 1` identifies the current
equation and discretisation contract. Analysis requires this explicit tag;
untagged data and earlier development versions are unsupported. Conservative
transport belongs to `polarisation_drift` version 1, selected explicitly with
`implementation=conservative`. Its `polarisation_ac_boundary`,
`polarisation_corrections` and `laplacian` settings belong to the drift section.
Completed parallel current enters auxiliary/vorticity evolution, and parallel
electrical work uses sheath-conditioned fields. The slab budget coverage remains
restricted; the version tag does not establish tokamak/sheath power-ledger
coverage. Reynolds and chi wall faces are retained for both AC conditions.
Numerical diagnostic names and conversions are preserved. Restart wall-storage
versioning is separate.

The auxiliary solve excludes the viscous rate, so polarisation transport
carries ion viscous drift work. `DivEviscous` remains diagnostic stress work;
it is excluded from the energy RHS and outward energy budget to avoid double
counting. Electron heating remains an internal transfer. Saved chi takes
precedence over reconstruction. The default is ten defect corrections.

Current continuity is `DivJpol+DivJdia+DivJpar+DivJextra+DivJcol=Rchi`, including
wall cells, with optional currents included when enabled. The existing X-Z
analysis paths are tested; general nonzero-parallel-current RHS/chi
reconstruction remains unfinished. The residual controls the gauge dependence
of local energy partition. Global budgets alone do not establish radial boundary
truncation accuracy. The exact wall rate would include `dt(psi_w(z))`. Chi retains
the electrostatic DC wall rate and default Neumann AC, neglecting this local
pressure-potential rate. The optional Dirichlet AC condition sets a zero AC
midpoint; it does not estimate the missing pressure-potential rate.

Use `PYTHONPATH=tools python -m energy_analysis.spectra RUN_DIR --output bands.json`
for fixed-mode wall/interior RMS, radial profiles and e-folding widths. This
streaming monitor supports singleton-Y, X-decomposed native version-1 dumps;
supply physical domain lengths if they differ from 50 mm. Bands 16–31 and
32–63 are fixed physical modes, independent of resolution. Wall amplitudes use
first physical cell centres. Widths include their half-cell offset from the
wall; crossings without wall dominance are flagged. A missing crossing has
no width. Native pytest regressions run with `HERMES_EXECUTABLE=build/hermes-3`;
add `--run-slow` for the drift-plane envelope test.

Optional numerical parallel phi damping uses collective sound speed. The
completed vorticity rate includes `ion_phi_dissipation`; auxiliary chi continuity
excludes both damping and viscosity. The diagnostic is a vorticity rate with the
`ddt(Vort)` conversion, not an external power sink.
The resolved option is saved as `ion_phi_dissipation_enabled`.

The wall condition uses conducting, Z-uniform phi, Neumann pressure and psi
guards, local Phi wall work and local species particle, momentum and
thermal/parallel kinetic weights. The native frozen-state and varying-fraction
regressions are in `tests/python/test_local_wall_transport.py`; nonlinear
convergence and regenerated wall envelopes use the fixed-wall D = 0.1 scan input.
