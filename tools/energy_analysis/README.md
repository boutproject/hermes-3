# Energy inventory and closed-slab balance analysis

Independent inventory auditing for `evolve_vorticity_ion_energy`, multiple ions
and a common ion temperature. Calculations use the C++ finite-volume face and
corner formulas, including `g13`; outputs are labelled xarray datasets.

Use the workspace `.venv` with editable xhermes:

```sh
uv pip install --python .venv/bin/python -e external/xhermes -r tools/energy_analysis/requirements.txt
.venv/bin/python tests/integrated/2D-ion-energy-curvature/analysis.py RUN_DIRECTORY --json energy-report.json
.venv/bin/python -m pytest -q tests/python external/xhermes/tests/test_energy_loading.py
.venv/bin/python tests/integrated/2D-ion-energy-equilibrium/runtest --exe build/hermes-3
.venv/bin/python tests/integrated/2D-ion-energy-curvature/runtest --exe build/hermes-3
.venv/bin/python tests/integrated/2D-ion-energy-multispecies/runtest --exe build/hermes-3
```

The runner creates a new temporary directory and prints its location. It keeps
inputs, simulation logs, exit status and JSON reports. Use `--output-root` to
select an existing parent directory. Source-tree and build-tree execution are
supported; CTest uses `HERMES_ENERGY_PYTHON` (default `.venv/bin/python`) and
registers each ion-energy directory as a separate CTest, plus
`energy-analysis-python`. MPI must be able to initialise its local transport
even for the single-process equilibrium test.

For Python use, add `tools` to the import path and call:

```python
from energy_analysis import open_run, energy_inventory, compare_diagnostics

run = open_run("RUN_DIRECTORY")
inventory = energy_inventory(run)  # joules, physical cells, time in seconds
comparisons = compare_diagnostics(run)  # absence of a diagnostic is explicit
```

`open_run` also accepts an xhermes Dataset. It reads evaluated options through
`bd.options`, retains physical guards for reconstruction, and excludes them
from integration. Both normalised and explicitly annotated SI inputs are
supported without modifying the caller's dataset. Unknown units are errors.
Path loading chunks one time record and Y plane at a time by default; pass
`chunks` to `open_run` to override this or supply an already chunked dataset.
Species names and charges come from resolved options; a `Species` manifest can
override discovery, while model and boundary options remain required.

The required saved fields are `Eion`, `Tion`, `phi`, `Vort`, `Pe`, charged species
densities, flowing ion momenta or prescribed velocities, and geometry
`J, dx, dy, dz, Bxy, g11, g33, g13`. Electron kinetic energy contributes to the
plasma inventory when `evolve_momentum` includes electron inertia. Algebraic
`zero_current` closure does not introduce an evolved electron kinetic reservoir.

The two ion inventories are saved `Eion` and reconstructed thermal plus parallel
plus perpendicular energies. They are compared without adding ion kinetic
energy to `Eion` twice. Optional diagnostics are read only by
`compare_diagnostics`, and never required by the inventory path.

Stage one supports slab topology, identity parallel transformation, named
unsqueezed `x,y,z` dimensions (including singleton Y/Z), and densities above
20 times both species and total ion recovery floors. Saved momentum is divided
by the selected component's soft-floored density. Cases near floors, missing
radial guards or unsupported tokamak topology fail explicitly. Perpendicular
operators support nonuniform and Z-dependent metrics and lazy/chunked arrays.
Physical integration assumes the model's `rho_s0**3 * J*dx*dy*dz` volume;
no full-torus multiplier is applied.

`energy_balance` checks that saved times are finite, increasing and regularly
spaced, then uses cumulative trapezoidal integration. It returns cumulative and
interval residuals, the reference energy and supplied power provenance. The
first interval residual is undefined (NaN, represented as null in JSON).
It requires explicit source and outward powers in watts for a complete residual;
`reservoir` selects an inventory field and defaults to `plasma_saved`. General
CLI reports only inventories and an incomplete balance until these
contributions are supplied. The equilibrium test uses analytically zero powers
and verifies its uniform, source-free, no-flow assumptions.

`electron_curvature_work` reconstructs `-Pe*Curlb_B.z*DDZ(phi)` from snapshots,
using scalar analytic curvature, saved normalisations and uniform saved `dz`.
It supports explicit, unfiltered FFT DDZ; it rejects grid-file curvature and
other schemes. The FFT derivative includes BOUT's zero Nyquist response and
rechunks complete Z lines while retaining lazy time/Y calculations. Use
`integrated_power` to convert a normalised energy RHS density to watts.
`polarisation_energy_coefficient` independently reconstructs multispecies
`E_adv`, including kinetic energy and the actual mass-density denominator;
its flux reconstruction is still pending.

## Integrated tests

Each total-ion-energy test owns a `data/BOUT.inp`, `README.md` and executable
`runtest`. The driven tests additionally own `criteria.json` with fixed
acceptance limits, cadence settings and source/boundary provenance. Their
shared audit code is in `tests/integrated/ion_energy_balance.py`.

| Test and CTest name | Coverage |
| --- | --- |
| [2D-ion-energy-equilibrium](../../tests/integrated/2D-ion-energy-equilibrium/README.md) | Uniform two-ion inventory, zero flow energy, diagnostics on/off |
| [2D-ion-energy-curvature](../../tests/integrated/2D-ion-energy-curvature/README.md) | Closed curvature exchange, cumulative/interval balances, cadence refinement |
| [2D-ion-energy-multispecies](../../tests/integrated/2D-ion-energy-multispecies/README.md) | Unequal masses/charges, varying composition, parallel energy and E_adv |

The individual READMEs explain each physical setup and power balance, exercised
terms, fixed tolerances, reports and limitations. The two driven `runtest`
scripts default to the `ci` profile: three actual cadence runs with diagnostics
on and an off run at the finest cadence. `--profile smoke` uses two cadences
and the off run. Both compare histories and strip optional outputs to verify
independent inventory and exchange reconstruction. Each test creates isolated
results with input hashes, run logs, exit status and machine-readable assertions.

`cadence_comparison` requires nested saved times at the same duration and
estimates both fine trapezoidal sampling error and the remaining residual by
second-order extrapolation. A third cadence checks the smooth-power assumption.
General source, open-boundary and flux reconstruction remain in the
[delivery plan](../../tests/integrated/ion-energy-analysis-plan.md).
The original [2D-energy](../../tests/integrated/2D-energy/README.md) remains the
separate `vorticity`/`evolve_pressure` interchange test.

Reports include per-record partition errors, diagnostic availability, forward
vorticity residuals, interpreter/import paths and revision provenance. Initial
records are not discarded. The driven tests exposed a restart-hook bug: the
first fresh-run output restored an unsolved initial `phi` after potential
recovery. Restart restoration now runs only at the first hook call, independently
for each component instance. C++ regressions cover fresh startup and multiple
instances, and the tightly solved integrated cases pass all records. Potential
recovery still has a finite iteration tolerance, so diagnostic comparisons use
explicit tolerances rather than assuming exact equality.
