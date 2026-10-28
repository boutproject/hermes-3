# Independent ion-energy inventory and power balance analysis

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
`reservoir` selects an inventory field and defaults to `plasma_saved`. The default CLI reports inventories and an incomplete balance. `--fluxes`
adds supported snapshot fluxes, and `--sources FILE` supplies an explicit SI
source manifest for complete supported balances. Missing sources or unsupported
geometry leave coverage explicitly incomplete. The equilibrium test uses analytically zero powers
and verifies its uniform, source-free, no-flow assumptions.

`electron_curvature_work` reconstructs `-Pe*Curlb_B.z*DDZ(phi)` from snapshots,
using scalar analytic curvature, saved normalisations and uniform saved `dz`.
It supports explicit, unfiltered FFT DDZ; it rejects grid-file curvature and
other schemes. The FFT derivative includes BOUT's zero Nyquist response and
rechunks complete Z lines while retaining lazy time/Y calculations. Use
`integrated_power` to convert a normalised energy RHS density to watts.
`polarisation_energy_coefficient` independently reconstructs multispecies
`E_adv`, including kinetic energy and the actual mass-density denominator;
`reconstruct_fluxes` applies this coefficient to the common polarisation flux.

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
Broader geometry, open electron momentum, sheath, conduction and reaction
coverage remain in the [delivery plan](../../tests/integrated/ion-energy-analysis-plan.md).
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

## Stage-three flux and source audits

Five independent directories provide the next budget tests:

| Test and CTest name | Physical coverage |
| --- | --- |
| [2D-ion-energy-fluxes](../../tests/integrated/2D-ion-energy-fluxes/README.md) | Driven unequal-ion ExB/Reynolds/electric/polarisation fluxes, sparse solve, species partition, regional Gauss budgets |
| [2D-ion-energy-heating](../../tests/integrated/2D-ion-energy-heating/README.md) | Known electron pressure source with thermal 3/2 conversion |
| [2D-ion-energy-sources](../../tests/integrated/2D-ion-energy-sources/README.md) | Analytic particle injection/removal, momentum forcing and thermal deposition |
| [2D-ion-energy-parallel](../../tests/integrated/2D-ion-energy-parallel/README.md) | Evolved electron inertia, pressure work, parallel MC transport and a quantified numerical-heating defect |
| [2D-ion-energy-boundary](../../tests/integrated/2D-ion-energy-boundary/README.md) | Open rectangular Y, imposed ion velocity, enthalpy/kinetic and electric boundary powers |

The new runners share `ion_energy_flux_test.py` and use the same input/log/report
contract as the closed tests. Supported operators live in `transport.py` and
`fluxes.py`. The current complete flux path requires saved uniform orthogonal
metrics (`g13=g23=g_23=0`), simple periodic or all-open rectangular Y topology,
identity transform, unfiltered FFT Z, C2 parallel work derivatives and MC
parallel limiting. It rejects unsupported components and dissipation explicitly.
This restriction does not narrow the nonuniform/nonorthogonal inventory and
Kperp operator support. Flux reconstruction additionally needs saved `g_22`,
`g23`, `g_23`, and any flowing electron closure momentum/velocity.

```python
from energy_analysis import (
    open_run,
    PrescribedSources,
    reconstruct_fluxes,
    reconstruct_sources,
    compare_flux_diagnostics,
    flux_power,
)

run = open_run("RUN_DIRECTORY")
# Fill from known external rates; each species, including e, needs an entry.
zero = {s.name: 0 for s in run.species}
sources = PrescribedSources(zero, zero, zero, "validated source-free input")
fluxes = reconstruct_fluxes(run, sources=sources)
checks = compare_flux_diagnostics(run, fluxes)
ion_parallel_outward_W = flux_power(run, fluxes, "DivEparallel")
```

`PrescribedSources` carries thermal [W/m³], density [m⁻³/s] and momentum
[N/m³] rates plus provenance. Scalar rates or labelled fields on normalised
run coordinates are accepted; all species require explicit entries. Total ion
source power includes thermal plus `v*Sm - A*v²*Sn/2`. Polarisation sources are
internal transport and are excluded. The source fixture uses a small component
compiled only with `HERMES_TESTS`, exposing scalar SI settings and no extra
production arrays. General spatial sources cannot be inferred from resolved
`<Field3D>` options; provide their known definitions rather than assuming zero.

Returned `*_flux_x/y/z` fields are **already integrated over a cell face**, in
normalised units, oriented towards increasing coordinate and indexed by the
lower-side cell. `flux_power` converts net outward faces to W, including explicit
contiguous subvolume slices such as `selection={"x": slice(2, 10)}`. No extra
area factor is applied. Energy faces are labelled `quantity=energy`;
`net_face_flux` returns normalised net currents or species fluxes, which
`flux_power` rejects to prevent labelling a current as watts. The spectral electrical Z divergence has an exact face
antiderivative with an arbitrary periodic constant-flux gauge; selected-volume
balances are independent of that gauge. Outer guard results are undefined.

The polarisation inverse is independently assembled as a sparse finite-volume
matrix. It uses a gauge for periodic X, or homogeneous half-face Dirichlet for
nonperiodic X. With `phi_boundary_relax=true`, the auxiliary χ retains
homogeneous midpoint Dirichlet guards; fixed-boundary runs retain the legacy
post-solve Neumann overwrite. It is
lazy across time/Y, with complete local X-Z solve planes. Diagnostic comparison
removes the free inversion gauge and reports both field and integrated errors.
The MC operators weight the common flux after limiting, so species charge
fluxes sum even with spatially varying composition.

`boundary_relaxation(run)` reconstructs radial potential midpoints, adjacent
Z-averaged targets and normalised rates from saved phi guards and `bd.options`.
Core Y averages exclude guards. The separately reconstructed
`DivEpol_enthalpy`, `DivEpol_parallel_kinetic` and `DivEpol_electric` face families
expose each energy contribution's physical-face closure. Particle and momentum
polarisation faces use the same mask. Full-current inversion permits normal
current; closed transport deliberately omits it at physical boundaries.

The diagnostic named `ddt(Phi)` contains an auxiliary current potential solving
`L(chi)=physical ddt(Vort)`, rather than the exact derivative of `phi+psi`.
Saved-output finite differences of `phi+psi` should not be compared as exact
χ diagnostics. Ordinary output `Vort` contains physical vorticity. New restart
files carry a version marker: their solver `Vort` physical radial boundary
cells encode normalised phi midpoints, which must not receive vorticity unit
conversion. Legacy restarts initialise those midpoints once from saved phi.

The separate `2D-ion-energy-relax-boundary` test checks diagnostics on/off and
stripping, output refinement, X/Y MPI decompositions and restart behaviour.
`tokamak-ion-energy-relax-boundary` uses six COMPASS ranks and a bounded raw-rank
audit of inventories, relaxation, sheath sources, parallel boundary power and
topology interfaces. It explicitly reports incomplete general tokamak power
coverage; the public slab flux reconstruction still rejects shifted topology.

For compact CLI reports:

```sh
PYTHONPATH=tools .venv/bin/python -m energy_analysis RUN_DIRECTORY --fluxes --json flux-report.json
PYTHONPATH=tools .venv/bin/python -m energy_analysis RUN_DIRECTORY --sources prescribed-sources.json --json balance-report.json
```

The source JSON has `thermal`, `density`, `momentum` dictionaries keyed by every
species and a `provenance` string, matching `PrescribedSources`. `flux_report`
keeps scalar ion/plasma histories, Gauss errors and diagnostic comparisons.
Unsupported electron kinetic ExB/open-boundary contributions leave the plasma
balance incomplete. A complete balance describes term coverage, not a claim
that its residual is zero.

## Tracked formulation/discretisation findings

- XPPM closes the periodic radial seam when `bndry_flux=false`; current
  reconstruction matches it. Global conservation does not validate seam fidelity.
- Reynolds current and energy fluxes use `Div_a_Grad_perp_XZ`, consistently
  with potential recovery, kinetic energy and polarisation transport.
  Reconstruction includes Z fluxes independently of the Y guard width.
- The parallel reconstruction and electron regression still track the former
  BOUT numerical viscous heating defect: twice the required wave-speed
  coefficient relative to the split particle/momentum flux. The electron
  case retains its raw nonconserving plasma residual and independently predicts
  the excess power with `numerical_heating_defect`. A separate `1e-8` relative
  defect-accounted gate validates the reconstruction, while the report explicitly
  marks raw plasma conservation false. BOUT's operator is now corrected;
  this reconstruction and regression still need updating.

The last finding means the electron parallel case is a regression/audit of a
known defect. Its README explains the coefficient identity and the meaning of
its pass criteria; it is not a source-free conservation gate.

## Pressure adapter and publication studies

`open_run` also accepts a slab using `vorticity` and saved species pressures.
`EnergyRun.formulation` identifies `pressure` or `total_energy`. Pressure runs
expose `ion_native` and `plasma_native` without inventing a saved `Eion`.
`ion_perpendicular` and `plasma_reconstructed` use a common ion-only pressure
potential and the FV estimator. Pressure native energy uses its actual signed
species pressure potential and the FD-X/spectral-Z quadratic form, including
radial boundary terms and the Z Nyquist mode. The pressure adapter currently
requires uniform orthogonal metrics for this native estimator.

`pressure.legacy_viscous_work` reconstructs the legacy corner-gradient heating
and its native kinetic work, exposing the heating mismatch. A corrected pressure
residual remains unavailable until legacy polarisation exchanges and full
`Phi` boundary-state work are covered. `reconstruct_fluxes` now includes scalar
kinematic viscosity for the total-energy slab, its face-weighted energy flux
and the compensating ion/electron heating transfer.

The [publication runner](../../tests/integrated/2D-energy/publication-study.md)
provides matched startup, streaming inventories, configurable resolution and
boundary matrices, conservative coordinate-aware restriction, PDF figures,
scalar tables and temporal/cadence pilots.

`closed_slab_outward_power` evaluates global outward power in the explicitly
closed singleton-Y, zero-parallel-flow contract without solving for interior
polarisation transport. It reconstructs Reynolds/viscous energy faces and uses
the implemented closure of the other normal faces; locally nonzero transport
remains available through the full reconstruction.
