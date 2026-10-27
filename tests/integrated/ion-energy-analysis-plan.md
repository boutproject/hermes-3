# Integrated tests and reusable energy analysis

Develop `2D-ion-energy-curvature/analysis.py` into a thin entry point for a reusable Python package,
proposed at `tools/energy_analysis/`. Use xhermes for loading and xarray for
labelled calculations. Start with the total ion energy model, common ion
temperature and slab geometry. Keep the pressure-based `data` case as a separate
comparison; it uses a different model and discretisation.
Use the workspace `.venv` and its editable `external/xhermes` installation for
all Python inspection, implementation, analysis and testing.

The first deliverable is an energy inventory and a cumulative balance that work
with `diagnose=false`. Subsequent deliverables reconstruct fluxes and sources,
compare them with optional diagnostics, and provide a more detailed account of
energy transfer. Stages one through three are now implemented, including the inventory
package, operator tests, diagnostic pairs, closed conservation tests and
independent slab flux/source/boundary audits; see [the module README](../../tools/energy_analysis/README.md). Conduction, reaction, broader geometry and tokamak budgets remain planned.

## Starting points verified in this checkout

- `2D-energy/runtest` runs the pressure-model interchange input `data`.
  Its pressure-sign assertions do not establish a quantitative energy balance.
  The total-ion-energy equilibrium, curvature and multispecies tests now have
  separate `2D-ion-energy-*` directories, inputs, READMEs and `runtest` scripts.
  The older `data-Eion` input is retained as `2D-ion-energy-curvature/reference`.
- The original `analysis.py` opened `data-Eion`, printed volume integrals and required
  `ddt(Eion)` and `ddt(Pe)`. It therefore fails when those diagnostics are absent.
- The previous short runtime checks covered slab, electron momentum, ion
  momentum, and two ions with different mass-to-charge ratios. Their inputs
  can seed checked-in cases. They used one output interval of 0.01 normalised
  time, so they verified execution rather than meaningful energy conservation.
- `Eion`, `Tion`, `phi` and `Vort` are written independently of the ion energy
  component's `diagnose` option. Evolved densities and momenta are also saved.
  `quasineutral` writes its density and `zero_current` writes its velocity and
  momentum independently of `diagnose`. Individual ion pressures are generally
  absent in the total ion energy model.
- xhermes opens physical units by default and reads `BOUT.settings`. Explicit
  `unnormalise=False` preserves the fields needed by the C++ discretisation.
  The resulting dataset exposes the parsed options through `bd.options`,
  including boolean conversion and expression evaluation methods.
- Before stage one, xhermes inferred species from short pressure-field names.
  On the two-ion output inspected here it reports only electrons. Species
  discovery must use resolved component settings and charges, not that inference.
- CMake currently copies only `runtest` and `data/BOUT.inp` for an integrated
  test. Shared analysis code and additional case inputs need explicit handling.
- `.venv/bin/python` resolves xhermes to `external/xhermes/xhermes/__init__.py`.
  pytest is installed in this environment for the stage-one tests.

## Data and units contract

Support a path or an already loaded xarray Dataset. The loader retains physical
guard cells for stencil calculations, and the integration mask excludes them.
It uses explicit spatial dimensions and mesh topology, rather than treating
every dimension other than time as a spatial dimension.

| Required information | Source and use |
| --- | --- |
| Time, `Eion`, `Pe` | Aggregate ion and electron thermal inventories; energy history |
| `Tion`, each ion density | Recover individual pressures with the common temperature |
| `phi`, `Vort` | Independent perpendicular energy and forward-operator check |
| Each flowing species' momentum or velocity and density | Parallel kinetic energy; include electrons when they have inertia |
| `J`, `dx`, `dy`, `dz`, `g11`, `g33`, `g13`, `Bxy` | Volume weights and the X-Z face/corner operators |
| Species names, `AA`, charge, active components | Species manifest, current closure and source coverage |
| Normalisations, density floors, solver and boundary options | Reproduce model conventions and convert results to physical units |
| Grid/topology and restart provenance | Physical boundaries, periodicity and unique time intervals |

Use `bd.options` as the options interface for a dataset `bd`; xhermes already
handles reading the options. Use section methods such as `get_bool` for boolean
options, `evaluate_scalar` for floating point expressions, and `evaluate` where
spatial expressions are needed. Do not add a separate settings parser or
expression evaluator. Read effective options, including saved defaults, through
this interface, and allow an explicit species manifest when options are
unavailable. Recognise arbitrary species names; punctuation and name length
must not determine charge or mass. Missing ions must produce a data-contract
error.

Use normalised fields for scheme reconstruction. Present inventories in joules
and powers in watts, using saved conversion factors and normalisation metadata.
For the present normalisation, physical cell volume is
`rho_s0**3 * J * dx * dy * dz`, and the energy density scale is
`Pnorm = e * Nnorm * Tnorm`. Check this against xhermes' physical geometry
conversions. Do not apply normalisations twice or insert a full-torus factor
into a simulated slab or toroidal wedge implicitly.

An adapter should also accept explicitly physical datasets by producing a
normalised calculation view. Test equivalence between the two routes. Do not
guess the units of a partially converted dataset from variable names.

In the common-temperature model, reconstruct normalised `p_s = n_s * Tion`.
For velocity and parallel kinetic energy, distinguish a saved physical species
momentum from an unmodified solver momentum. Audit dump timing and the
`evolve_momentum` density-floor correction before promising equivalence near
the floor. Reproduce `softFloor(n, f) = max(n, 0) + f*exp(-max(n, 0)/f)` wherever
the selected component uses it; a hard floor is not interchangeable. Initially
use densities well above floors for conservation gates, then add floor cases
with a separately validated reconstruction contract.

Keep required halo values. A missing radial guard cell cannot generally be
replaced with a periodic roll. Reconstruct guards only when the actual boundary
rule and its parameters are known. Otherwise report a supported interior
subvolume or an unavailable whole-domain reconstruction. Retain singleton Y
dimensions and handle both `MYG=0` and nonzero guard widths.

## Independent inventory and operator reconstruction

Implement Python operators from the documented face/corner formulas, without
calling the C++ implementation. Use named dimensions, periodic Z shifts and
radial halos. Avoid a cell-centred `numpy.gradient` approximation for the
scheme-matched perpendicular energy.

In normalised variables, reconstruct

\[
\psi = \sum_s \frac{A_s p_s}{\bar A Z_s},\qquad
\Phi = \phi + \psi,\qquad c = \frac{\bar A}{B^2}.
\]

Reimplement the X-Z forward operator and `calculateKperp` using the existing
diagonal face fluxes and symmetric `g13` corner terms. For integrated face
fluxes `F`, the cell energy is

\[
K_{\perp,i} = \frac{1}{4J_i dx_i dz_i}
\sum_{\text{four faces}} \Delta\Phi_{\rm face} F_{\rm face}.
\]

Return the face fluxes as well as the divergence and energy. This permits a
direct integration-by-parts check with radial boundary work, and avoids deriving
boundary power by subtracting unrelated large volume integrals. The full cell
volume also includes `dy`; the X-Z flux formulas are for each Y plane.

Expose two ion inventories:

1. The saved total `integral(Eion)`.
2. The independently reconstructed `integral(3/2*sum(p_s) + sum(Kparallel_s)
   + Kperp)`.

Compare these to test the pressure/kinetic partition. Add `3/2*Pe` and, where
applicable, electron parallel kinetic energy to obtain the plasma inventory.
Do not add ion kinetic energies to `Eion` again. A pressure-based model gets its
own adapter and operator selection; its spectral-Z inversion must not be
declared equivalent to this finite-volume discretisation.

Validate operators on synthetic xarray datasets before relying on simulations:
constant fields, Fourier modes, nonuniform metrics, nonzero `g13`, radial
boundary work, periodicity, variable coefficients and operator symmetry.
Cross-check small simulation outputs against `ion_KE_perp`, `ion_KE_parallel`
and `ion_pressure_potential` when available. Use tight recovery tolerances in
these cases, since the diagnostic kinetic energy is evaluated during the
pressure/potential iteration, not necessarily at precisely the final iterate.

## Balance and flux analysis

Define positive source power as entering the plasma and positive boundary power
as leaving it. Make the primary history check

\[
R(t) = E(t)-E(t_0)-\int_{t_0}^{t}Q_{\rm ext}\,dt'
                         +\int_{t_0}^{t}P_{\rm out}\,dt'.
\]

This works without saved RHS diagnostics. Report absolute and normalised
residuals, the energy scale used, cumulative source/boundary contributions and
coverage. Use interval balances as well as the final value so opposing errors
cannot hide through cancellation.

Time integration assumes regular outputs. Before integrating, check that the
time axis is finite, strictly increasing, has at least two samples, and has
constant spacing within a documented floating point tolerance. Derive the
spacing from the saved time values rather than the configured output interval.
Reject irregular spacing, gaps or duplicate restart records with a clear error;
do not silently integrate them using a single timestep. Use cumulative
trapezoidal integration and estimate its uncertainty by output-cadence
refinement. Snapshot differences are interval rates, not the solver's
instantaneous RHS; use centred derivative estimates only with adequate samples,
and avoid endpoint derivative estimates as tight conservation gates.

Develop the flux reconstruction in this order:

| Terms | Reconstruction and validation |
| --- | --- |
| X-Z perpendicular operator | Shared diagonal/cross face geometry, volume divergence and boundary work |
| ExB and Reynolds transport | Reproduce the actual stream-function interpolation, limiting, sign conventions and boundary switches |
| Electrical work | Ion `phi*Ji`, electron electrical energy flux and thermal ExB compression; track species exchanges separately |
| Parallel transport | Species enthalpy and kinetic fluxes, face geometry, numerical wave speeds and boundary-flux rules |
| Polarisation transport | Common current, linearly interpolated species weights, actual mass-density denominator and energy carried per current |
| Sources and conduction | Prescribed sources first; then supported heat exchange, friction, conduction and reactions |

Reconstruct the model's normalised polarisation energy coefficient

\[
E_{\rm adv} = \phi + \sum_s\frac{A_s}{Z_s\rho}
                  \left(\frac52p_s + K_{\parallel s}\right),\qquad
\rho=\sum_s A_s n_s,
\]

including the configured mass-density floor. Apply weights to the common
reconstructed Reynolds face flux, rather than limiting every species' weighted
field independently. Check that the species charge fluxes sum to the common
current above the floor.

For scheme-level polarisation cross-checks, implement an independent sparse
finite-volume solve for `ddt(Phi)` from the reconstructed implemented vorticity
RHS. Match the solve boundary flags, gauge and subsequent transport boundary
rules, including explicit closure of radial `g13` fluxes. A temporal derivative
of reconstructed `Phi` is also useful when RHS reconstruction is unavailable,
but label it as a cadence-dependent estimate; boundary handling can make it
different from the transport potential used by the code.

Return integrated face powers with their orientation and units documented.
Do not multiply an already area-integrated flux by area a second time. Check
discrete Gauss balance for each supported flux family and cancellation across
periodic and processor interfaces. Provide boundary maps and selected
subvolume budgets on request.

With diagnostics enabled, compare reconstructed fields and integrals with the
saved `DivE*`, `DivJpol`, `phi_Div_Je`, `electron_exb_work`, `Eion_src` and time
derivatives. Also check the instantaneous ion RHS identity. Electron thermal
and kinetic budgets remain separate: `phi*Div(Je)` is not an electron thermal
source. Include the electron momentum work before assessing total electron or
plasma energy.

External total-energy sources include thermal energy plus
`v_s*S_momentum_s - 1/2*A_s*v_s**2*S_density_s`. Internal polarisation transport
already belongs to the flux and must not be counted again as an external
source. Track collision exchanges between species so internal transfers do not
appear as net plasma heating.

State snapshots cannot recover every arbitrary reaction source, numerical
flux or boundary history exactly. Each term must carry provenance such as
`reconstructed`, `diagnostic`, `prescribed`, `estimated` or `unavailable`.
Never replace an unavailable source or boundary term with zero. Offer a complete
balance when all needed contributions are known, and otherwise report an
explicitly incomplete balance and unexplained power.

Audit the existing fields and `evolve_vorticity_ion_energy::outputVars` before
adding outputs. Reconstruct inventories, fluxes and prescribed sources from
already available fields and `bd.options` wherever possible. If an essential
contribution cannot be recovered, the default outputs may be extended outside
the `diagnose` block, keeping additions to the minimum required for power
balance auditing. Prefer static metadata or compact scalar integrated sources
and boundary powers over new time-dependent volume fields. Record the reason
each addition is necessary and its storage cost. Keep detailed spatial energy
flow diagnostics optional.

Validate any new summaries against independent Python calculations and
diagnostic runs before using them as references. A summary from this component
must state whether it covers only ion energy or the full plasma budget; electron
and other component contributions need their own accounting. Closed/source-free
and prescribed-source tests should work without additional budget outputs.

## Integrated case matrix

Promote the previous runtime inputs into small checked-in slab cases. Start
with identity parallel transformation and fixed potential boundaries; add more
geometry and physics only after the simpler balances can be interpreted.

| Case | Main assertion |
| --- | --- |
| Uniform equilibrium | Constant inventory, common ion temperature, zero flow energies and zero fluxes |
| Curvature-driven X-Z slab | Measurable electron/ion energy transfer; independent inventories and cumulative balance including any boundary power |
| Electron parallel momentum | Include electron kinetic energy and electrical work; avoid a false thermal-only conservation check |
| Ion parallel momentum | Species enthalpy/kinetic transport and particle/momentum source conversion |
| Two ions with unequal `A/Z` | Common temperature, spatially varying composition, current partition and multispecies `E_adv`; reference mass density different from actual density |
| Nonuniform, nonorthogonal slab | Radially varying widths, nonzero `g13`, forward residual and integration by parts with actual saved metrics |
| Prescribed thermal heating | Inventory increase equals the known deposited energy; begin with the existing electron pressure source and its 3/2 conversion |
| Particle injection/removal | Known thermal, momentum and particle source energy accounting; use a test-only prescribed-source component if no existing component supplies the controlled combination |
| Boundary transport | Independently measure outward powers, including Reynolds/electrical work, rather than assuming zero flux from a potential boundary label |
| COMPASS tokamak, six MPI ranks | Global and regional power balances with core, wall and divertor boundary fluxes; topology and processor-interface cancellation |

For parallel cases, use explicitly periodic Y for the first closed-domain
checks. Add open ends and sheath/conduction cases separately with their flux
accounting. For nonorthogonal geometry, provide a consistent metric tensor and
disable orthogonal metric recalculation; verify the saved `g13`, `J` and `Bxy`
actually match the intended geometry.

### Tokamak boundary stress test

Use the checked-in `examples/tokamak-2D/compass-36x48.grd.nc` for a small
axisymmetric total ion energy case (`MZ=1`), with electrons evolved using
`evolve_pressure` and a common ion temperature. Use the existing tokamak
examples as geometry and boundary setup references, replacing their
pressure-based ion model with `evolve_vorticity_ion_energy`. Run on six MPI
ranks, allowing oversubscription; for Open MPI the launcher is
`mpirun --oversubscribe -np 6`. Validate a decomposition compatible with the
grid topology and record it in the case manifest. Stage the grid in the
isolated run directory so the test does not rely on the working directory or
the examples' `tokamak.nc` convention.

Load the actual saved metrics, topology, guard cells and parallel-transform
settings. Build physical boundary masks for the core, outer wall, private-flux
radial boundaries and divertor targets as present in this grid. Periodic core
connections, separatrix branch connections and MPI interfaces are internal
faces, and their contributions must cancel rather than count as boundary power.
Check that each physical cell and boundary face is included exactly once.

Start with controlled initial profiles, supported parallel transport and
known core input or prescribed heating, with measurable power leaving through
open targets. Add sheath heat losses and conduction once their discrete flux
accounting is available; add recycling and reactions later. Keep each case's
source and boundary coverage explicit. Report a whole-plasma balance and
core/SOL/private-flux subvolume balances where the topology supports them,
including cancellation of transfers between regions. Include electrical and
polarisation boundary work as well as thermal and parallel kinetic transport.

Run diagnostic on/off pairs, cadence refinement and discrete Gauss checks.
Assert nonzero input and target losses, independent-versus-diagnostic agreement
for supported terms, and bounded cumulative residuals after sampling and solver
errors are resolved. Compare alternative compatible decompositions when
practical. Add this case to a bounded integration profile after the slab flux
tests pass; include a six-rank oversubscribed smoke run with a runtime limit.

### Run criteria and refinement

Run conservation candidates long enough to produce a measurable transfer or
transport signal, not just startup output. Require the exercised contribution
to be comfortably larger than the allowed balance error. The earlier one-step
runs remain useful smoke checks but must not be the conservation gate.

Every core case has paired runs with all relevant component diagnostics enabled
and disabled. Compare physical histories as well as analysis results. Also
remove optional diagnostic variables from an enabled dataset in a Python test
to prove the inventory path does not read them accidentally. Those two checks
test different failure modes and both are needed.

Use independent parameter sweeps rather than an unrestricted Cartesian matrix:

- Output cadence at fixed mesh, physical duration and solver tolerances.
- Spatial resolution at fixed physical domain, wavenumbers and duration.
- Potential recovery and time-solver tolerances independently.
- One and two MPI ranks, including a split that exercises processor face fluxes.
  The tokamak test additionally uses six oversubscribed ranks and a
  topology-compatible decomposition.
- Nontrivial `Nnorm`, `Tnorm`, `Bnorm` and explicit normalised/physical input
  views; reconstruct the same physical energy through either view.
- Standard metric builds first, followed by a metric-3D case with Z-dependent
  geometry and then shifted-coordinate tests.

Known model or boundary inconsistencies should become named tracked findings,
with their measured contribution and refinement behaviour. Do not loosen a
single total-energy tolerance until they disappear from the report, or assert
that all current terms are already implemented because a startup case passes.

## Module interface and test integration

### Development environment and xhermes changes

Use `.venv/bin/python` directly for standalone analysis and
`.venv/bin/python -m pytest` for Python tests. Configure CMake/CTest to use this
same interpreter. Keep the editable xhermes installation in place so changes in
`external/xhermes` are exercised immediately. Verify the interpreter and
xhermes import path in test setup and record the Hermes and xhermes revisions
in reports, together with whether either checkout has local changes.

Inspect and amend the xhermes source as needed for the analysis. Initial areas
to review are:

- `xhermes/load.py`: species and charge discovery from `bd.options` or explicit
  field metadata, including total ion energy runs without individual ion
  pressure outputs; conversion attributes and dimension assumptions.
- `xhermes/accessors.py`: normalised-versus-physical geometry, singleton
  dimension handling, guard widths and topology metadata. In particular,
  geometry extraction currently squeezes dimensions and gives derived `dv`
  and `da` physical-unit labels; validate those labels when
  `unnormalise=False`.
- `xhermes/selectors.py` and existing selector/decomposition tests: core, SOL,
  private-flux and target selections for the COMPASS case, with and without
  guards and across processor layouts.

Put reusable loader, metadata, geometry and diagnostic corrections in xhermes;
keep formulation-specific finite-volume operators and energy budgets in the
energy-analysis package. Add regression tests alongside any xhermes changes,
including pressure-based loading, total-energy multispecies loading, both unit
representations and guard handling as applicable. Reuse the checked-in COMPASS
grid for offline tests: the existing xhermes geometry tests download other
fixtures, which should not become a prerequisite for the local core suite.
Select additional diagnostic changes from concrete reconstruction or
cross-check failures and apply the minimal-output policy above.

### Package and runner

Proposed package responsibilities:

| File | Responsibility |
| --- | --- |
| `dataset.py` | xhermes loading, `bd.options` access, normalisations, species manifest and domain/guard contract |
| `operators.py` | Pure face/corner operators, kinetic energy and discrete identities |
| `inventory.py` | Ion/electron energy densities and volume inventories |
| `fluxes.py` | Supported flux families, current split and source reconstruction |
| `balance.py` | Cumulative/instantaneous balances, coverage and uncertainty |
| `report.py` and `__main__.py` | CLI, compact JSON/NetCDF summaries and optional plots |

Representative public functions are `open_run`, `calculate_kperp`,
`energy_inventory`, `reconstruct_fluxes`, `energy_balance` and
`compare_diagnostics`. Return labelled datasets with units and provenance,
not printed arrays or figures as the only result. The curvature test's `analysis.py`
is a small wrapper that accepts a path and report options; importing the
package performs no file loading or analysis.

Keep xarray calculations lazy where possible. Stream time records and Y planes
for X-Z operators; require complete local halos or implement overlap explicitly
for spatially chunked data. Cross-chunk results must match unchunked results.
Reduce to scalar time series by default and retain full reconstructed fields or
boundary maps only when requested. Static geometry should not be duplicated
at every time step. Import plotting libraries only for plotting.

Replace the current runner's pressure-sign criterion for the total-energy
variant with assertions from a machine-readable report. Preserve the legacy
pressure test under its existing model identity. Resolve the executable through
the build/test configuration, accept an executable override, and create a new
isolated output directory for each case. Record exit status and logs before
analysis; do not reuse stale output or delete user simulation data.

Extend the CMake copy/install logic for the package, case manifests and COMPASS
grid fixture. Register the tokamak case with six MPI processes, explicit
oversubscription support and a runtime limit. Exercise both source-tree and
build-tree test execution, including paths containing spaces. Declare the
xhermes, xarray, NumPy, SciPy and pytest dependencies explicitly; add missing
dependencies to the existing workspace `.venv` and retain the editable xhermes
checkout. Use this environment for CTest and standalone
analysis. Keep generated results and fixtures separate from checked-in inputs.

## Delivery order and acceptance

1. **Inventory foundation.** Verify the workspace `.venv` and editable xhermes,
   add missing test dependencies, and make required xhermes metadata/geometry
   corrections with regression coverage. Implement the loader and species/units
   contract, Python perpendicular operators and ion/electron inventories. Add synthetic
   operator/unit tests and replay small existing outputs. Acceptance: analytic
   identities pass, enabled diagnostics agree within recovery/roundoff error,
   options expressions are handled through `bd.options`, and optional
   diagnostics can be removed without breaking inventory analysis.
2. **First conservation tests (implemented).** Check in equilibrium, curvature and two-ion
   cases; run actual diagnostic on/off pairs. Add cumulative balances and JSON
   assertions. Acceptance: inventories agree, exercised terms are nonzero where
   intended, boundary/source coverage is explicit, and cadence refinement
   separates sampling error from model residuals. Irregular or invalid time
   axes must fail the regular-output precondition.
   The separate curvature and multispecies `runtest` CI profiles run at
   intervals 12, 6 and 3 over duration 120, plus finest-cadence
   diagnostics-off pairs. Reconstructed FFT curvature work closes the electron
   budget without diagnostics; plasma drift is below `4e-16` of reference
   energy, and electron residuals fall by about four per cadence halving.
   Fixed tolerances and source/boundary provenance are in each driven test's
   `criteria.json`; the audit implementation is shared in `ion_energy_balance.py`.
   The restart-hook initial-phi discrepancy is fixed and covered by unit tests.
   The periodic radial ExB seam remains a tracked transport-fidelity finding;
   it does not break these global closed-domain balances.
3. **Flux and source budgets (implemented for uniform orthogonal slabs).** Add ExB/Reynolds/electrical, parallel and
   polarisation reconstruction in that order, followed by controlled sources
   and boundary cases. Acceptance: discrete Gauss checks, diagnostic field and
   integrated comparisons, and species current summation all pass for supported
   schemes. Missing terms produce an incomplete report rather than a false pass.
   Separate `2D-ion-energy-fluxes`, `-heating`, `-sources`, `-parallel` and
   `-boundary` directories now own input, README, criteria and runners. The
   module reconstructs MC ExB/parallel transport, Reynolds/electrical and
   weighted polarisation fluxes, an independent sparse FV `ddt(Phi)` inverse,
   and prescribed thermal/particle/momentum source conversion. Diagnostic
   field/integral, ion RHS identity, species charge face/current, Gauss and
   two-subvolume checks run with enabled/off/stripped diagnostics. The open
   boundary input is a simple rectangle with imposed ion velocity, not a sheath.
   No production outputs were added; the controlled-source component is test-only
   and retains scalar SI rates in saved options.

   Flux coverage is explicitly limited to uniform orthogonal identity slabs,
   unfiltered FFT Z, C2 parallel work and the saved MC limiter. Unsupported
   sources remain unavailable without an explicit provenance-bearing manifest.
   Independent Kperp/forward operators retain their nonuniform/g13 coverage.

   Three findings remain named: XPPM's closed periodic radial seam, BOUT Reynolds
   FV omission of Z when MYG=0, and excess parallel numerical heating from its
   wave-speed coefficient. The electron parallel test matches diagnostics and
   independently predicts the last defect. Its report preserves the raw plasma
   drift and marks source-free conservation false; a separate defect-accounted
   residual gate checks reconstruction without loosening the conservation limit.
   Open electron momentum, sheath/conduction/reaction operators and general
   kinetic ExB fluxes need additional discrete contracts in subsequent stages.
4. **Broader integration.** Add the six-rank COMPASS tokamak boundary test,
   MPI and resolution/tolerance sweeps, metric-3D and shifted geometry. Extend
   default audit outputs only where reconstruction has an identified gap.
   Acceptance: tokamak global and regional balances account for physical
   boundary powers, internal interfaces cancel, and budgets and coverage are
   reproducible across layouts, unit representations and diagnostic settings;
   remaining physical discrepancies are quantified.

Use a small smoke profile, a bounded CI profile for the core matrix, and a full
profile for refinement and geometry studies. Choose fixed tolerances from the
initial convergence studies, then check them into the manifests. Do not choose
them independently for each run. Suggested starting targets for tightly solved
small cases are normalised absolute operator errors around `1e-12`, relative
diagnostic inventory differences around `1e-9`, and cumulative balance errors
around `1e-6` of the reference energy. These are provisional, not claims about
the present full model: accept them only after resolving solver, quadrature,
boundary and known discretisation contributions.

The resulting reports should answer: how much energy is stored in each
reservoir, what entered or left the selected volume, what moved between species,
how much remains unexplained, and which information or approximation limits
that conclusion.
