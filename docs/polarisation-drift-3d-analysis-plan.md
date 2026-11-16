# Polarisation drift and 3D ion-energy analysis: agent implementation plan

## Objective and agreed scope

Move the conservative polarisation drift calculation out of
`evolve_vorticity_ion_energy` and into `polarisation_drift`, then validate total
ion energy evolution with parallel currents and sheath boundaries. Deliver an
xhermes-based analysis of a 3D tokamak turbulence example that accounts for
input power, physical X and Y boundary power, and domain energy.

When assigned this plan, the implementation agent should carry that first
task through the checkpoints below. The agent must write a review report on
completion or on the first failed checkpoint. This planning turn adds only
this document; it does not implement the code changes.

Required outcomes:

1. An enum selects `legacy` or `conservative` polarisation drift.
2. Existing X-Z ion-energy behaviour is preserved, including conservative
   currents, source conversion, wall transport and energy budgets.
3. A 3D MMS variant exercises nonzero parallel current and its divergence.
4. A production-sheath regression verifies boundary/current ordering.
5. A separate 3D tokamak example, comparable to
   `examples/tokamak-3D/tcv-x21`, evolves ion total energy.
6. Its analysis reports thermal, parallel kinetic and perpendicular kinetic
   energy, actual input power, physical boundary powers and a power-balance
   residual, with explicit coverage and numerical uncertainty.

Defer Reynolds extraction into separate files, all changes to `vorticity.cxx`
and `vorticity.hxx`, and the conservative pressure-versus-total-energy
comparison. Conservative drift initially pairs with
`evolve_vorticity_ion_energy`; legacy drift pairs with unchanged `vorticity`.
Keep Reynolds evaluation in the ion-energy component. Do not change the
Reynolds discretisation to accommodate the refactor.

## Review reports and checkpoint policy

### Maintainer follow-up: phi damping, wall convergence and chi acceptance

The [ion-energy damping and convergence plan](ion-energy-dissipation-convergence-plan.md)
governs the newly authorised follow-up. Add optional `phi_dissipation` only
to physical `ddt(Vort)` in ion-energy `finally()`, leaving chi, polarisation
transport and the energy RHS unchanged. Investigate the Reynolds wall failure
with the current MC-plus-R form and remove obsolete production-attribution
diagnostics such as bracket Q. Do not change production Reynolds/wall methods
or the legacy `vorticity` component in this follow-up.

The new parallel-current MMS's absolute `ion_chi <= 0.001` cap becomes
diagnostic-only; retain chi errors and spatial/temporal/MPI convergence checks.
This supersedes references below to all six stage-3 gates remaining unchanged
only for that absolute chi assertion. Other thresholds and named deferrals
remain in force. Preserve the original raw six-failure evidence and report
re-evaluation under the new policy separately. The follow-up plan specifies
the bounded targeted diagnostics allowed for the known convergence failures;
unexpected implementation failures and failed coupled validation still
require a report and stop of dependent work. Stages 4–5 remain pending.

### Maintainer amendment after stage-1 review

Checkpoint 0 is accepted, with the parallel kinetic exchange minimum revised
to `5e-5`; energy/recovery and other exchange thresholds remain unchanged.
The existing `test_reynolds_wall_curvature_converges` accuracy failure is
explicitly deferred for this task. Keep the test enabled, its `>1.5` Reynolds
order threshold unchanged, and record its measured failure separately at each
relevant regression check. This named existing limitation does not stop the
remaining preservation, parallel-current/sheath, MMS or tokamak work. Any new
failure or worsening of the preserved Reynolds behaviour still follows the
stop/report policy below. Do not change the Reynolds method, extract Reynolds
into separate files, or change `vorticity.cxx`/`vorticity.hxx` to address it.

The stage-1 constructor-order defect must be fixed and covered with drift
listed before ion energy, both with default and explicit drift inverse types.
Complete the interrupted wall and unstarted MMS preservation checks before
declaring checkpoint 1 passed and proceeding to the parallel-current/sheath
work. This amendment defers the named Reynolds test and the eight existing default
MMS accuracy gates listed below; all other acceptance limits and required
outcomes remain in force.

### Maintainer amendment: existing default MMS accuracy gates

Accept the constructor fix. Defer only these eight `ion-energy-mms` default
8-to-16 order gates for this task:

- `order/baseline/Ni`, `order/baseline/Pe`, `order/baseline/Vort`,
  `order/baseline/Eion`;
- `order/curvature/Ni`, `order/curvature/Pe`, `order/curvature/Vort`,
  `order/curvature/Eion`.

Keep each test enabled and its `>=1.7` threshold unchanged; retain raw failures
and classify this named deferral separately from passing assertions. All eight
historical `criteria.json["observed"]` measurements already fell below 1.7.
Six fresh pre-extraction results reproduce them closely, but baseline and
curvature vorticity decline from about 1.69 to 1.38. Investigate that decline
separately; this amendment does not attribute it to Reynolds or authorise a
numerical-method change. The older passing-suite narrative is inconsistent
with the recorded criteria. These failures precede the extraction.

Complete all remaining MMS and wall preservation checks before proceeding.
Bitwise preservation and all other budget, recovery, field-error, diagnostic
and restart gates remain required. This deferral excludes long-profile gates
and the new parallel-current MMS, which must pass its own accuracy requirements.

### Maintainer amendment: observed stage-2 ordering roundoff

Accept the retained production-sheath reversed-construction trajectory
comparison as a harmless ordering effect: maximum normalised absolute
physical-field difference `8.881784197001252e-16`. The logged schedules differ
in ion-momentum placement but retain recovery -> sheath -> conservative drift;
identical-state C++ RHS and initial integrated fields agree bitwise.
Retain the 14 raw `sorting/*` failures and their zero thresholds in the
original evidence. This acceptance covers those observed results only; it does
not relax new MMS accuracy, physical current/electrical checks, MPI limits or
same-input X-Z preservation. Finish stage-2 X-Z preservation, then continue the
plan if all remaining checkpoint requirements pass. The precise source of the
ordering roundoff remains unisolated; no permission dependency or numerical
method should be changed merely to remove it.

Write the report to `docs/polarisation-drift-3d-analysis-report.md`. Keep logs,
measurement JSON, plots and run directories under a documented validation
directory, preferably `build/validation/polarisation-drift-3d/`. Link retained
artifacts from the report. Do not overwrite unrelated user files or earlier
run evidence.

At each checkpoint, record the commands, input/build identifiers, assertions,
thresholds and measured results. On the first failed checkpoint:

- Write the report immediately with status `CHECKPOINT FAILED`, naming the
  failed stage and assertion, expected result, measured result, and logs.
- Stop dependent stages and hand the report back for review.
- Preserve the failing inputs and outputs. Explain any established cause and
  the smallest proposed next step; distinguish evidence from hypotheses.
- Do not relax thresholds, suppress failing terms, declare missing powers
  zero, change boundary physics, or perform an unplanned numerical-method
  change to make the check pass.

An unavailable required grid, dependency, diagnostic or reference dataset is
a failed prerequisite check, not a completed validation. Report the missing
item and exactly which outcomes could not be tested. Routine local edits and
build preparation are allowed before a checkpoint is evaluated; once the
checkpoint fails, report and hand off.

On completion, write status `COMPLETE` only when all required outcomes have
passed. The report must contain:

- Scope and changed files, including any deviations from this plan.
- Hermes, BOUT++ and xhermes revisions; executable/library hashes; resolved
  input settings and grid identity; Python environment and MPI decomposition.
- A checkpoint table with measurements and acceptance limits.
- The 2D before/after comparison and legacy regression results.
- MMS orders and nonzero-term checks; sheath ordering/current checks.
- Tokamak run duration, resolution, enabled physics and demonstrated
  turbulence, including whether a stationary regime was reached.
- The power ledger, energy reservoirs, instantaneous and cumulative
  residuals, coverage, and timestep/output-cadence sensitivity.
- Remaining limitations and proposed follow-up work. Do not present a short
  startup run as a validated long-duration turbulence simulation.

Report a newly encountered baseline failure too. Identify known historical
failures separately: `docs/test-consolidation-validation.md` records an
optional long MMS skew-vorticity order failure. Historical measurements do
not substitute for checks of the present checkout, and this task must not
weaken that existing threshold.

## Source map and initial audit

Read the applicable repository instructions and inspect the working tree
before editing. Preserve unrelated changes, including `external/BOUT-dev`
and the existing local `external/xhermes` checkout.

Primary implementation locations:

- `src/evolve_vorticity_ion_energy.cxx`, especially recovery, the polarisation
  block beginning near line 1287, finalisation and source conversion.
- `include/evolve_vorticity_ion_energy.hxx`.
- `src/polarisation_drift.cxx`, `include/polarisation_drift.hxx`.
- `src/sheath_boundary.cxx`, relevant parallel transport/source components.
- `tests/unit/test_evolve_vorticity_ion_energy.cxx` and
  `tests/unit/test_polarisation_drift.cxx`.
- `tests/integrated/ion-energy-mms`, `ion-energy-closed`,
  `ion-energy-sources` and `ion-energy-walls`.
- `examples/tokamak-3D/tcv-x21/BOUT.inp`.
- `tools/energy_analysis`, `tests/python`, and the publication scan's existing
  tests in `examples/tokamak-2D-driftplane/2D-energy`.

Known starting constraints to verify rather than assume:

- The current MMS varies in X, Y and Z but uses `zero_current` for electrons;
  the total parallel current is cancelled.
- Ion-energy finalisation currently omits parallel current divergence from
  the vorticity rate, while some parallel electrical terms are evaluated
  before sheath conditioning.
- Conservative transport stores private polarisation source records to avoid
  counting their kinetic conversion twice.
- The analysis loader currently rejects tokamak topology and non-identity
  parallel transforms; sheath losses and core-only sources have incomplete
  budget coverage.
- TCV-X21 evolves electron density and derives ion density through
  quasineutrality. Its ion pressure component owns the configured ion heat
  source, which disappears when that pressure component is removed.

Check the prior parallel-work and sheath follow-ups in
`docs/ion-energy-follow-ups.md`. Where historical implementation is useful,
inspect it through git. Reuse xhermes topology support rather than restoring
a second independent rank-file/topology reader.

## Stage 0: establish baselines and the data contract

Before changing equations, capture representative current X-Z RHS/diagnostic
outputs and the existing integrated energy checks. Include multiple species,
nonuniform/cross metrics, both auxiliary AC boundary conditions, fixed and
relaxing radial walls, and viscosity on/off using existing coverage. Preserve
a legacy drift regression with unchanged `vorticity`.

Audit which numerical fluxes, powers and inventories are already saved and
which are needed for the 3D budget. Write a field contract listing each
variable's owner, physical meaning, units/conversion, dimensions, face index
and orientation, boundary ownership, and whether it is required for analysis.
Separate implementation diagnostics from mandatory budget outputs.

Specify acceptance thresholds before comparing changed results. Same-input
2D refactoring should agree to roundoff where arithmetic is unchanged; any
tolerance for reordered operations must be justified by baseline precision,
not fitted to a regression. Preserve existing MMS and integrated thresholds.

**Checkpoint 0:** Baseline checks pass and the data contract can account for
all terms in the proposed model. Report a missing prerequisite or unexpected
baseline failure before proceeding.

## Stage 1: move conservative polarisation drift

Add a BOUT-style enum `implementation = legacy | conservative`. Default to
legacy to preserve existing standalone drift configurations; set conservative
explicitly in every migrated ion-energy configuration. Legacy keeps the
current potential-flow, compression and option behaviour. It does not need
legacy Reynolds face reconstruction in this task.

The ion-energy component must publish the inputs needed by conservative
drift: Boussinesq coefficient, pressure potential, Reynolds face currents and
their divergence, diamagnetic-current divergence, and auxiliary DC wall
rates. These inputs must be independent of `diagnose`. Keep the Reynolds
calculation and its numerical method in the existing component.

Move the auxiliary potential solver, fixed-count defect corrections,
polarisation currents, weighted species transport and wall averaging into
`polarisation_drift`. Move auxiliary-specific types/helpers as needed, but do
not use this as an opportunity to factor or redesign Reynolds stress.

Conservative drift must publish the energy-flux divergence and separately
identified particle/momentum source contributions consumed by ion-energy
finalisation. Replace the private cross-component source cache with state
fields. Preserve the existing total-energy source conversion and the exclusion
of polarisation's already-transported parallel kinetic energy.

Preserve:

- `Jpol = -JReynolds - c grad_perp(chi)` with the current code's face/sign
  convention and its existing poloidal extension.
- Actual species mass-density partition, density floors, pressure-potential
  weighting and uniform wall products.
- The fixed number of spectral solves, zero start on each RHS evaluation,
  homogeneous correction boundaries and physical wall currents.
- The exclusion of the kinematic-viscous rate from the auxiliary RHS and the
  resulting viscous drift work carried by polarisation transport.
- Deterministic RHS evaluation and accepted-output recovery warm starts.

Migrate affected inputs and analysis field ownership/version handling. Preserve
useful diagnostic names or document their mapping; diagnostic ownership may
change, physical 2D behaviour must not. Reject unsupported component pairing
with a clear configuration error. Do not add conservative thermal-pressure
accounting or alter `vorticity` in this task.

**Checkpoint 1:** Legacy regressions pass. Current and migrated X-Z
ion-energy runs agree within the predeclared limits for recovery, derivatives,
particle/momentum sources, polarisation currents, energy transport, wall
exchanges and integrated budgets. Run the appropriate affected C++ and Python
tests and existing default energy profiles. Document any changed metadata.

## Stage 2: parallel currents and production sheath ordering

Use component permissions to establish:

`recovery -> sheath boundary -> conservative drift -> all finally() methods`.

Read charged-species momentum with its boundary region valid. Form parallel
current from `sum_s (Z_s/A_s) momentum_s`, retaining momentum rather than
reconstructing `N*V`. Include its divergence in both the auxiliary RHS and
the ion-energy vorticity rate. Treat enabled extra/collisional current fields
consistently between those equations; retain the existing viscous convention.

Defer parallel electrical flux/work to completed-state reads, including the
sheath-conditioned potential. Verify parallel slices, communication,
field-aligned transforms and the existing poloidal transport extension.
Do not require sheath-conditioned momentum in the recovery stage, which would
reintroduce the dependency cycle.

For the chosen current conventions, check the identity

`DivJpol + DivJdia + DivJpar + DivJextra + DivJcol = Rchi`,

where `Rchi` is the actual selected auxiliary RHS minus its forward operator;
omit absent optional fields, not enabled physical terms. Check physical
boundary cells as well as the interior. Periodic perpendicular planes require
a solvable mean RHS; do not silently discard a nonzero mean current source.

Add a focused open-Y regression with the actual sheath component, independently
check the imposed current at each target, and verify that changing momentum
guards changes the current divergence and auxiliary result as predicted.
Test automatic component sorting and an MPI decomposition exercising Y
communication. Keep MMS boundary prescriptions separate from production
sheath conditions unless the manufactured solution satisfies them exactly.

**Checkpoint 2:** Sorting is acyclic; sheath current, electrical work and
auxiliary continuity checks pass; the X-Z baseline remains unchanged.

## Stage 3: 3D parallel-current MMS

After the checkpoint-3 accuracy review, the maintainer authorised a separate
electron-momentum RHS audit and tighter 16-cell temporal comparison before
prioritising the 32-cell diagnostic run. Keep these diagnostic continuations
separate from the default runner, which must evaluate each completed
refinement pair and stop through its result-writing context before launching
subsequent validation on failure. The six new accuracy gates and all other
thresholds remain unchanged. Stages 4–5 remain pending.

Extend the existing symbolic continuum MMS machinery with a parallel-current
variant. Manufacture electron momentum independently of ion momentum, with
nonzero total current and nonzero divergence. Add the evolved electron
momentum reference/forcing as appropriate; do not retain `zero_current` in
this variant.

Derive forcing from the intended continuum equations. Construct an analytic
auxiliary inverse including parallel current that satisfies the chosen X
conditions; start with periodic Y to isolate equation correctness. Ensure
that Reynolds, parallel-current and polarisation terms have finite measurable
magnitudes and are not engineered to cancel every term being tested.
Check MMS forcing ownership so it is injected once, not also interpreted as
a physical heating source.

Use at least three spatial refinement levels with sufficiently tight temporal
and potential-recovery errors. Measure L2 and Linfinity errors for evolved
fields and recovered pressure, temperature, phi, psi and chi. Check expected
spatial orders and the auxiliary residual; repeat one level with tighter
temporal settings to establish separation of errors. Include an MPI check.
Register a bounded default test and document any longer profile.

**Checkpoint 3:** Symbolic identities and nonzero-term checks pass, followed
by the specified convergence and temporal-sensitivity checks. Save exact
fields, forcing definitions, criteria and measured orders in the report.

## Stage 4: xhermes-based analysis implementation

### Reuse and upstream-compatible structure

Use xhermes features whenever possible. The scripts remain here for now, but
their reusable functions are intended to become part of xhermes.

Use `xhermes.open` / `open_hermesdataset` for loading output and resolved
settings, with the grid supplied for toroidal geometry. Reuse `.hermes`
accessors, SI conversion, species metadata, radial/poloidal region selectors,
topology/target metadata and geometry support where applicable. Inspect the
installed API rather than assuming every local helper supports 3D. The local
checkout has several geometry routines with explicit X-theta array shapes;
verify their handling of the additional toroidal dimension before using them.

Core budget functions should accept xhermes-loaded `xarray.Dataset` objects
and return `Dataset`/`DataArray` results with units, provenance, coverage and
missing-term metadata. Keep plotting, file paths, study execution and CLI
configuration outside the numerical core. Existing `EnergyRun` APIs may remain
as thin adapters for current tests; do not create a new parallel data model
for tokamak analysis.

Support labelled dimensions, including xhermes toroidal/poloidal naming,
broadcasting and chunked arrays. Avoid positional reshaping or eager loading
of full time histories. Compute reductions explicitly when writing reports.
Do not mutate the caller's dataset during unit conversion. Use package imports,
not a hard-coded repository path to `external/xhermes`.

Add only missing functionality locally. Isolate minimal 3D adapters/extensions
and document their proposed upstream API. Do not edit the external xhermes
checkout as part of this task. Where an upstream feature is inadequate, retain
explicit coverage checks and report the gap rather than copying an independent
loader, species parser or topology implementation.

### Inventory, domain and boundary contract

Integrate only physical cells using the grid Jacobian and cell volumes. Use
xhermes region/topology information to exclude guard cells and branch-cut
duplicates and identify physical X walls and Y targets. Periodic/inter-rank
connections are internal faces, not losses. Validate these masks against the
solver's physical boundary ownership.

Distinguish a toroidal cell volume, a simulated-sector volume and a full-torus
quantity. Apply any sector-to-torus factor consistently to every source,
boundary power and inventory. Report both the factor and domain convention.
Do not apply an extra 2*pi factor to already integrated face powers or volumes.

Report in SI units:

- Ion and electron thermal energy: `integral 3/2 P_s dV` for this model.
- Ion and electron parallel kinetic energy: `integral 1/2 m_s n_s V_s^2 dV`,
  using the solver's recovered velocity/floor convention.
- Ion perpendicular kinetic energy from the same face operator as recovery,
  counted once for the Boussinesq ion flow.
- Saved Eion and its reconstruction from the three ion reservoirs.
- Total domain energy `Eion + Ue + Kparallel_e`, compared with the independent
  reservoir sum; do not add ion reservoirs to Eion a second time.

Prefer saved solver-consistent perpendicular energy for the tokamak budget
when independent reconstruction would require an unsupported shifted-grid
operator. Label its provenance and retain an independent reconstruction check
on supported test geometry. Do not approximate tokamak perpendicular energy
with an orthogonal Cartesian gradient.

### Power ledger and required numerical diagnostics

Implement a ledger before summing a residual. For every enabled energy term,
classify it as external input, physical boundary exchange, internal transfer
or an explicit model/numerical defect. Identify its owner and reservoir.
Internal electron-ion exchanges must cancel in the plasma sum; do not assume
mechanical/electrical work cancels when an evolved reservoir is absent.

Use the actual numerical fluxes used in the RHS. Add the minimum solver output
support required by the ledger, without changing its discretisation. Small
boundary/integral histories may be sufficient; expose detailed face diagnostics
where needed to verify ownership and signs. Every reported scalar must have a
documented integration/normalisation and a check against its numerical term.

Required entries include:

- Electron/ion heating, source masks and time prefactors actually applied.
- Particle/momentum injection or removal energy conversion, distinguished
  from internal polarisation transport.
- Signed inner/outer X boundary powers: enabled ExB, diamagnetic/electrical,
  polarisation, diffusion/conduction and other energy-carrying terms.
- Signed power for each physical Y target: numerical parallel enthalpy and
  kinetic transport, electrical work and additional sheath heat corrections.
- Any configured radiation, explicit sinks, wall-potential work or numerical
  energy defects needed by the selected model.

Sheath heat corrections are currently stored as adjacent-cell energy sources.
Identify their contribution separately and count them as boundary power,
not external heating. Conduction may close boundary flux because sheath heat
is already imposed; account for what the solver actually does. Do not add a
gamma-based total sheath power on top of advective energy and its additional
correction. Likewise, do not add diagnostic `DivEviscous` as an extra energy
loss when polarisation transport already carries that drift work. Include
electrical and wall work exactly once.

An enabled term without a measured or justified reconstruction makes coverage
incomplete. Retain known histories, list missing terms and omit a claimed
closed residual. Output cadence uncertainty, RHS integration error and
auxiliary-current residual are separate quantities, not interchangeable
explanations for missing physics.

### Budget results

With positive outward powers, report

`rP = dE/dt - Pin + PXout + PYout + Pother_loss`

and

`RE(t) = E(t) - E(0) - integral_0^t (Pin - PXout - PYout - Pother_loss) dt`.

Treat a diagnosed energy defect separately and show whether including it
explains the residual; do not hide it in Pin. Signed entry through a boundary
is allowed. Show power in watts, energy/residual in joules, and residuals
normalised by a declared input-power/integrated-throughput scale. Avoid a
near-zero instantaneous denominator.

Use the cumulative balance as the primary transient test. Where available,
also compare the spatial integral of the instantaneous RHS with the power
ledger, separating temporal quadrature error from spatial/source closure.
Report reservoir histories, per-boundary powers, internal-transfer checks,
recovery consistency and auxiliary-current residuals. Do not require
Pin = Pout during startup when domain energy is changing.

**Checkpoint 4:** Tests using xhermes-loaded or equivalently annotated
datasets verify SI/normalised equivalence, no input mutation, physical-domain
volume, target masks, face signs/ownership, sector scaling, source-versus-sheath
classification, missing-term behaviour and energy reconstruction. Exercise
multi-target topology and multi-rank collection. Preserve existing slab,
wall, source, native analysis and publication tests. Cross-check at least one
supported budget with solver face/integral diagnostics. Record the remaining
3D API gaps and how the local implementation handles them.

## Stage 5: tokamak example and physical power-balance validation

Create `examples/tokamak-3D/tcv-x21-ion-energy/` with documented input, grid
requirements, run instructions and a thin xhermes-based analysis entry point.
Use the same supplied tokamak grid as the reference case; fingerprint the
actual grid and record any resolution differences. If it is unavailable, do
not claim the tokamak validation passed.

Use deuterium, evolved ion density and parallel momentum, quasineutral
electron density, evolved electron pressure and parallel momentum,
`evolve_vorticity_ion_energy`, conservative drift and production sheaths.
Retain comparable shifted geometry, normalisations and input-power profiles.
Document why the density ownership differs from the existing TCV-X21 setup:
conservative ion particle transport must enter an evolved continuity equation.
Verify that dependent electron density, particle-source conversion and
parallel kinetic energy accounting remain consistent.

Provide explicit ion heating independent of ion pressure evolution, with SI
units and saved applied power. Match the reference's nominal species heating
where practical, but calculate actual applied power from the source field,
core mask and volume; do not trust a copied shape_factor. Construct initial
Eion and vorticity consistently from the intended physical fields.

Run a small 3D sheath smoke case first, then a finite-duration turbulence case.
Enable the reference's collisional, conductive and drift physics in controlled
increments, checking the ledger each time. Do not silently transplant
vorticity-only dissipation switches into the ion-energy component. Missing
support for an essential term requires a documented deviation or failure
report, not an untracked numerical replacement.

Compare at least a representative interval at tighter solver settings and
finer output cadence. Set defensible power-residual acceptance limits before
the final comparison using independent budget benchmarks and temporal
refinement estimates; do not fit limits to the tokamak residual. Include a
spatial/decomposition check of the terms most
sensitive to physical boundaries. Report the observed residual, not just a
pass flag. A residual comparable to input power, or one dominated by an
unidentified term, fails the budget validation.

Produce plots of all energy reservoirs, Pin and per-boundary Pout, transient
dE/dt and rP, cumulative RE, and the relevant recovery/current residuals.
Provide a machine-readable budget and coverage report. State whether the run
is transient or stationary and how turbulent fluctuations were established.

**Checkpoint 5:** The example advances with nonzero parallel current and
production sheath fluxes, exhibits documented 3D turbulence, has a complete
ledger for its enabled physics, and passes the predeclared power-balance and
temporal-sensitivity checks. Unsupported terms, unexplained deficits or a
missing required run are failures requiring review.

## Final validation and handoff

Use the repository's configured test environment; discover the actual build
directory and Python executable. Applicable existing checks include affected
unit fixtures, `energy-analysis-python`, the default `ion-energy-*` integrated
tests, and the publication scan's Python tests. Register and run the new
parallel-current/sheath checks. Avoid rerunning unrelated broad suites once
the required checks pass unless a change or unresolved failure justifies it.

Run formatting/lint and `git diff --check` appropriate to changed files. Ensure
that `vorticity` and Reynolds discretisation remain unchanged, and that no
unrelated external checkout changes were made. Update equation/configuration
documentation, the analysis README and deferred-follow-up notes to reflect
the implemented scope.

Write the completion report with all checkpoint evidence and links to the
new example, analysis API/CLI, numerical output contract and retained plots.
List future work separately: Reynolds extraction, selectable conservative
Reynolds in pressure evolution, a matched conservative thermal-pressure
equation, and upstreaming the reusable analysis functions into xhermes.
