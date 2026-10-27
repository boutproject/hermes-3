# Energy inventory and closed-slab balance analysis

Independent inventory auditing for `evolve_vorticity_ion_energy`, multiple ions
and a common ion temperature. Calculations use the C++ finite-volume face and
corner formulas, including `g13`; outputs are labelled xarray datasets.

Use the workspace `.venv` with editable xhermes:

```sh
uv pip install --python .venv/bin/python -e external/xhermes -r tools/energy_analysis/requirements.txt
.venv/bin/python tests/integrated/2D-energy/analysis.py RUN_DIRECTORY --json energy-report.json
.venv/bin/python -m pytest -q tests/python external/xhermes/tests/test_energy_loading.py
.venv/bin/python tests/integrated/2D-energy/runtest-inventory --exe build/hermes-3
.venv/bin/python tests/integrated/2D-energy/runtest-balance --exe build/hermes-3 --profile ci
```

The runner creates a new temporary directory and prints its location. It keeps
inputs, simulation logs, exit status and JSON reports. Use `--output-root` to
select an existing parent directory. Source-tree and build-tree execution are
supported; CTest uses `HERMES_ENERGY_PYTHON` (default `.venv/bin/python`) and
registers `2D-energy-inventory`, `2D-energy-balance` and `energy-analysis-python`. MPI must be able to
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
spaced, then uses cumulative trapezoidal integration. It returns cumulative and
interval residuals, the reference energy and supplied power provenance. The
first interval residual is undefined (NaN, represented as null in JSON).
It requires explicit source and outward powers in watts for a complete residual;
`reservoir` selects an inventory field and defaults to `plasma_saved`. General
CLI reports only inventories and an incomplete balance until these
contributions are supplied. The equilibrium test uses analytically zero powers
and verifies its uniform, source-free, no-flow assumptions.

Stage two adds `balance-curvature` and `balance-two-ion`, with fixed duration
120 in normalised time and actual runs at output intervals 12, 6 and 3. Both
have periodic X/Z, periodic core Y connections, identity transformation and
Y-uniform fields. The mixture has unequal masses/charges, varying composition
and opposite ion flows. Inputs prescribe zero pressure/density sources; a
component whitelist excludes unaccounted heating/momentum sources. Saved
topology, periodic X halos, Y uniformity and disabled dissipation are checked
before supplying zero external and outward powers. Detailed flux fidelity is
left to stage three.

`electron_curvature_work` reconstructs `-Pe*Curlb_B.z*DDZ(phi)` from snapshots,
using scalar analytic curvature, saved normalisations and uniform saved `dz`.
It supports explicit, unfiltered FFT DDZ; it rejects grid-file curvature and
other schemes. The FFT derivative includes BOUT's zero Nyquist response and
rechunks complete Z lines while retaining lazy time/Y calculations. Use
`integrated_power` to convert a normalised energy RHS density to watts.
`polarisation_energy_coefficient` independently reconstructs multispecies
`E_adv`, including kinetic energy and the actual mass-density denominator;
its flux reconstruction is still pending.

The runner checks electron thermal change against integrated reconstructed
curvature work, and plasma change against the known closed/source-free budget.
Enabled runs also compare `electron_exb_work` pointwise, `phi_Div_Je` after
integration, kinetic/pressure diagnostics, zero external sources and zero
integrated energy-flux divergences. The finest cadence has an actual diagnostics
on/off pair; their state and inventory histories must agree. Removing all
optional time-dependent outputs from the enabled dataset must preserve both
inventories and exchange reconstruction.

Each run retains `execution.json`, the input SHA-256, `run.log` and
`energy-report.json` with histories in joules, powers in watts, coverage and
named pass/fail assertions. The top-level `result.json` records cadence
comparisons and pair assertions. `--profile smoke` uses intervals 6 and 3 (six
runs); `ci` uses all three (eight runs). The separate equilibrium test remains
part of CTest. `cadence_comparison` requires nested saved times at the same
duration and estimates the finer trapezoid error; a third cadence checks its
smooth-power assumption.

The initial convergence study found electron residuals, divided by the plasma
reference energy, of `3.34e-9, 8.37e-10, 2.09e-10` for the curvature case and
`1.80e-9, 4.50e-10, 1.12e-10` for the mixture. The factor of about four per
halving identifies sampling error. Plasma drift was below `4e-16`, partition
error below `2e-15` normalised and forward vorticity error below `4e-11`.
The fixed manifest tolerances allow solver/roundoff margins; the minimum
exchange and flow signals are much larger than the balance limits.

Tracked finding **periodic radial ExB seam**: the current XPPM operator uses
`firstX/lastX` with `bndry_flux=false`, suppressing transport across the radial
periodic seam. It remains globally conservative in these cases. These tests
establish energy conservation and exchange; they do not establish correct
periodic-interface transport. Open boundaries, general source accounting,
parallel fluxes, MPI interface fidelity and tokamak auditing remain later work.

Reports include per-record partition errors, diagnostic availability, forward
vorticity residuals, interpreter/import paths and revision provenance. Initial
records are not discarded. The driven tests exposed a restart-hook bug: the
first fresh-run output restored an unsolved initial `phi` after potential
recovery. Restart restoration now runs only at the first hook call, independently
for each component instance. C++ regressions cover fresh startup and multiple
instances, and the tightly solved integrated cases pass all records. Potential
recovery still has a finite iteration tolerance, so diagnostic comparisons use
explicit tolerances rather than assuming exact equality.
