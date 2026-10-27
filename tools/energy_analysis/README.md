# Energy inventory analysis

Independent inventory auditing for `evolve_vorticity_ion_energy`, multiple ions
and a common ion temperature. Calculations use the C++ finite-volume face and
corner formulas, including `g13`; outputs are labelled xarray datasets.

Use the workspace `.venv` with editable xhermes:

```sh
uv pip install --python .venv/bin/python -e external/xhermes -r tools/energy_analysis/requirements.txt
.venv/bin/python tests/integrated/2D-energy/analysis.py RUN_DIRECTORY --json energy-report.json
.venv/bin/python -m pytest -q tests/python external/xhermes/tests/test_energy_loading.py
.venv/bin/python tests/integrated/2D-energy/runtest-inventory --exe build/hermes-3
```

The runner creates a new temporary directory and prints its location. It keeps
inputs, simulation logs, exit status and JSON reports. Use `--output-root` to
select an existing parent directory. Source-tree and build-tree execution are
supported; CTest uses `HERMES_ENERGY_PYTHON` (default `.venv/bin/python`) and
registers `2D-energy-inventory` and `energy-analysis-python`. MPI must be able to
initialise its local transport even for the single-process equilibrium test.

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
spaced, then uses cumulative trapezoidal integration. It requires explicit
source and outward powers in watts for a complete residual. General CLI reports
only inventories and an incomplete balance until these contributions are
supplied. The equilibrium test uses analytically zero powers and verifies its
uniform, source-free, no-flow assumptions. Flux reconstruction, driven
conservation tests and tokamak auditing belong to subsequent stages.

Reports include per-record partition errors, diagnostic availability, forward
vorticity residuals, interpreter/import paths and revision provenance. Initial
records are not discarded: pressure/potential recovery can leave a diagnostic
kinetic energy evaluated at the preceding iterate. On an older default-tolerance
multispecies fixture the initial perpendicular energy mismatch was about
`7.1e-5` normalised, while the next record agreed below `1e-14`. Treat that as a
recovery/timing discrepancy to investigate, rather than a Python conservation
pass. The tightly solved equilibrium passes all records.
