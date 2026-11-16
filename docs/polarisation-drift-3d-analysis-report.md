# Polarisation drift and 3D ion-energy analysis: review report

**Status: CHECKPOINT FAILED — stage 3 parallel-current MMS accuracy.**

Checkpoint 2 passed with the maintainer's acceptance of the observed sorting
roundoff. The completed 8/16-cell new MMS references now fail six predeclared
accuracy gates: vorticity L2 order, four Linfinity orders, and the chi L2
field-error cap. Of 66 evaluated checks, 60 pass. The existing MMS/Reynolds
deferrals do not apply. Thresholds remain unchanged; stages 4 and 5 have not
started. Fine-grid, temporal-sensitivity and MPI checks are incomplete.
The workflow correction now stops the default runner before subsequent
validation on failure. The separate tighter 16-cell diagnostic hit both its
600-second and extended 1800-second output watchdogs, so temporal sensitivity
remains unassessed.

The P2 constructor-order defect is fixed. All 80 affected C++ tests pass,
including reversed construction with default/explicit drift inverse types
and bitwise RHS agreement. The maintainer-authorised policy amendment defers
the known Reynolds wall-accuracy test while retaining its test and thresholds.
Checkpoint 0 remains accepted, with parallel kinetic exchange minimum `5e-5`.

Both default MMS profiles completed all five runs, with all 57 assertion values
identical: 49 pass and eight existing `>=1.7` order gates remain raw failures
under the narrow deferral. Both wall profiles completed all 13 executions and
passed all 204 assertions, with identical measured values. All 148 MMS and 343
wall physical field histories agree bitwise, including owned wall faces and
potential guards; 234 MMS and 566 wall numerical histories agree too. Timing
histories are excluded explicitly. Diagnostic-state differences are zero;
fixed/no-slip restart errors are `9.688780000516317e-9` and
`6.756289794296322e-8`, below `1.2e-6`.

Checkpoint 1 preservation has passed under the named existing-accuracy policy.
The following stage-1 evidence predates the stage-2 implementation and
failed checkpoint documented below. No Reynolds/vorticity method or MMS threshold has changed.

## Stage-2 review: observed sorting roundoff accepted

The completed production-sheath regression fails its predeclared bitwise
reversed-construction history gate: maximum absolute difference
`8.881784197001252e-16`, required `0`. All 14 compared physical histories differ
at later output times. This gate is not included in the existing MMS/Reynolds
deferrals. Checkpoint 2 originally failed and stopped dependent work. The
maintainer now accepts these observed differences as probably harmless
ordering effects. The raw 14 failures and zero thresholds remain retained;
this is an acceptance of those measurements, not a general accuracy relaxation.
Stage-2 X-Z before/after preservation has now passed: all 207 common numerical
histories agree bitwise, with every common assertion identical. Closed profiles
pass 41/41 current assertions; source profiles pass 15/15. The extra closed
assertion checks the newly saved `ddt(Vort)` diagnostic and passes at
`2.927345865710862e-17` versus `1e-10`. Stage-1 preservation
remains valid.

The identical-state C++ scheduler comparison passes bitwise. The three complete
integrated executions (ordered/reversed two-rank Y decomposition and ordered
serial) pass target-current, parallel-divergence, completed electrical-work
and auxiliary-continuity gates, including target-adjacent physical cells.
They also pass the serial/Y-MPI tolerance. Of 61 assertions, 47 pass and 14
history-sorting gates fail; raw results and thresholds remain unchanged.

| Stage-2 check | Acceptance | Measured result |
| --- | --- | --- |
| Affected C++ tests, including actual sheath, guard perturbation and periodic mean | All pass | 83/83 pass |
| Identical-state reversed scheduler RHS | Bitwise | Passed |
| Target current from independent Bohm/electron formulas | Absolute error `<=1e-12` | Maximum `3.552713678800501e-15`; passed |
| Current-divergence reconstruction | Absolute error `<=1e-12` | Maximum `0`; passed in all three runs |
| Auxiliary continuity / physical vorticity rate | Absolute error `<=1e-10` | Maximum `1.3877787807814457e-17` / `0`; passed |
| Completed ion/electron parallel electrical expressions | Absolute error `<=1e-11` | Maximum `0`; passed |
| Serial versus two-rank physical histories | Absolute error `<=1e-10` | Maximum `1.5543122344752192e-15`; passed |
| Reversed-construction physical histories | Bitwise, absolute difference `0` | **Failed:** maximum `8.881784197001252e-16` |
| Stage-2 X-Z before/after preservation | Existing bitwise requirement | Passed; 207 common histories bitwise, all common assertions identical |

The resumed [checkpoint-2 decision](../build/validation/polarisation-drift-3d/stage2-accepted-order-pf3ri46u/checkpoint-2.json)
and [X-Z preservation comparison](../build/validation/polarisation-drift-3d/stage2-accepted-order-pf3ri46u/preservation-comparison.json)
retain frozen inputs, before/current library and executable fingerprints,
commands and full closed/source results. The new diagnostic adds one assertion;
comparison matches common assertions by name, rather than assuming unchanged
list positions. No physical field differs. The initial comparator's mismatched
assertion-list positions were a validation-script issue, not a numerical
failure. Checkpoint 2 now passes with the maintainer's narrow ordering acceptance.

The retained `BOUT.log.0` execution lists **differ**: ordered runs execute
ion density -> electron quasineutrality/pressure/momentum -> ion momentum;
reversed runs execute ion density -> ion momentum -> electron
quasineutrality/pressure/momentum. Both then execute ion-energy recovery ->
sheath -> conservative drift. See the
[execution-order comparison](../build/validation/polarisation-drift-3d/stage2-r1fdwln9/component-execution-order.json).
The scheduler uses this same list for `finally()` methods. Ion momentum's
transform permissions read its own density and mass, so its placement before
electron components is permitted by the current declarations. A different
legal topological sort does not alone establish missing permissions, but
component/finalisation dependencies and shared-state effects must be audited
alongside solver registration order. The initial-zero-momentum unit comparison
does not exclude a dependency exposed by evolved states.

The failure does not establish an RHS ordering defect: the identical-state
production scheduler unit test passes, and all 14 integrated fields agree exactly
at the initial saved time. The [time-by-time comparison](../build/validation/polarisation-drift-3d/stage2-r1fdwln9/sorting-time-comparison.json)
records the differences arising during integration. Changed time-integrator variable
registration/reduction order is a hypothesis, not an established cause. The
smallest next step is to audit the differing component/finalisation order on
identical nonzero-momentum evolved states, then isolate solver registration
order if those RHS comparisons still agree bitwise. The maintainer review now accepts the observed roundoff; all other
checkpoint requirements and new accuracy gates remain required.

Stage 2 implements completed charged momentum reads, parallel-current
communication and divergence in the auxiliary/vorticity RHS, optional
extra/collisional currents, completed-potential electrical work, and a
volume-weighted periodic-plane solvability rejection without projection.
Ion-energy output version is now 6; conservative drift output version is 2.
The analysis recognises versions 4/5/6, including pressure-formulation paths.
Version-6 tag support is validated for existing X-Z paths only; general
nonzero-parallel-current RHS/chi reconstruction and a production-sheath power
ledger remain unfinished, so these runs do not claim full analysis coverage.
The production test is registered as `ion-energy-sheath` with two ranks.
Reynolds discretisation and `vorticity.cxx`/`vorticity.hxx` are unchanged.

Python verification: 112 analysis tests pass. The full Python run records
160 passes, one skip, the unchanged named Reynolds accuracy failure
(orders `0.0252637567`, `0.0048517793` versus `>1.5`), and three MPI startup
failures inside the sandbox. Retrying those native tests outside the sandbox
passes both gauge cases. The wall test initially encounters stale metadata
assertions for versions 5/1; updated assertions require implemented versions
6/2, and the wall-energy check then passes. No numerical acceptance limit is
changed. Thus 163 Python checks pass in aggregate, one is skipped and the named
Reynolds gate remains a raw failure. See
[full Python log](../build/validation/polarisation-drift-3d/stage2-r1fdwln9/python-all.log),
[native retry](../build/validation/polarisation-drift-3d/stage2-r1fdwln9/python-native-retry.log)
and [final native wall check](../build/validation/polarisation-drift-3d/stage2-r1fdwln9/python-native-version6.log).
Formatting, Ruff and `git diff --check` pass. No dependent stage was started.

The first integrated validation attempt used solver `NV` output as if it
retained target guards. In this BOUT revision, non-FCI `setBoundaryTo` skips
new Y boundary regions, so the restored solver output has zero initial Y
momentum guards. The first attempt also requested an unsaved vorticity-rate
field. That incomplete attempt is retained as test-preparation evidence,
not counted as a completed checkpoint validation. New `parallel_current`,
`ion_parallel_current`, `electron_parallel_current` and `ddt(Vort)` diagnostics
save the completed numerical values without changing momentum evolution.
The target check independently computes sheath formulas from physical cells;
the companion C++ test checks each species target momentum directly and the
analytic response to a guard-only perturbation. The periodic fixture was also
corrected to select identity transforms rather than the unrelated FCI mock.

Evidence directory:
[stage2-r1fdwln9](../build/validation/polarisation-drift-3d/stage2-r1fdwln9/).
See [checkpoint-2.json](../build/validation/polarisation-drift-3d/stage2-r1fdwln9/checkpoint-2.json),
[raw integrated result](../build/validation/polarisation-drift-3d/stage2-r1fdwln9/ion-energy-sheath-4echyuf9/result.json),
[sheath log](../build/validation/polarisation-drift-3d/stage2-r1fdwln9/sheath-final.log),
[C++ log](../build/validation/polarisation-drift-3d/stage2-r1fdwln9/unit-diagnostics.log)
and [C++ XML](../build/validation/polarisation-drift-3d/stage2-r1fdwln9/unit-diagnostics.xml).
The checkpoint JSON fingerprints the executable, library, unit binary and
new test input/runner. Every execution records its command, input hash,
executable hash and resolved settings; `NXPE=1`, `MYG=2`, `ny=8` exercises Y
communication on two ranks. All runs use `FI_SOCKETS_IFACE=en0` outside the
sandbox. The before-stage-2 source and executable/library are preserved.

```sh
cmake --build build --target hermes-3 hermes_unit_tests -j 4
FI_SOCKETS_IFACE=en0 build/hermes_unit_tests '--gtest_filter=*IonEnergy*:*PolarisationDrift*:*VorticityTest*'
FI_SOCKETS_IFACE=en0 .venv/bin/python tests/integrated/ion-energy-sheath/runtest --nproc 2 --output-root build/validation/polarisation-drift-3d/stage2-r1fdwln9
.venv/bin/python -m pytest -q tests/python/test_energy_analysis.py
```

## Stage-3 review: new parallel-current MMS accuracy failure

The new [parallel-current MMS](../tests/integrated/ion-energy-parallel-mms/README.md)
extends the existing continuum generator with independently evolved electron
momentum, removes `zero_current`, and adds a finite analytic parallel current
and auxiliary inverse. The original four generated inputs are byte-identical
under regeneration; all five variants pass symbolic preflight and file checks,
and all 11 symbolic/generated-input Python tests pass. The new test is
registered separately as `ion-energy-parallel-mms`; original default/long
selections and thresholds remain unchanged.

With `a=pi X/10`, `b=2 pi Y/10`, `theta=2 pi Z/10-t/5`, `tau=t/5`, `c=2`,
`Jpar=0.03 cos(2a) sin(b-tau) cos(theta)`, `Je=Jpar-Ji` and
`NVe=-Je/1836`. The added chi is
`-0.03 cos(2a) cos(b-tau) cos(theta)/(4 c pi/10)`.
Symbolic checks verify `c lap(chi)=DivJReynolds+DivJpar`, current formation from
momenta, nonzero-term magnitudes, periodic Y/Z and the Neumann X AC / zero DC
wall conditions. Full weighted polarisation divergences retain finite
`DivJpol=-DivJpar`; particle, momentum and energy forcing is derived from the
continuum equations, rather than a numerical reference. Electron momentum
forcing contains parallel inertia, pressure gradient, electrical force and
ExB advection. BOUT MMS adds it once; physical heating diagnostics are checked
separately.

Predeclared [stage-3 acceptance](../build/validation/polarisation-drift-3d/stage3-parallel-mms-4en7eogy/acceptance.json)
requires three spatial levels 8/16/32 (Y twice X/Z), L2 orders >=1.7 then >=1.8,
Linfinity orders >=1.5, auxiliary/continuity/recovery errors <=1e-10,
nonzero physical terms, a tighter 16-cell temporal run with differences <=5%
of spatial error (absolute floor 1e-11), and serial/Y-MPI agreement <=1e-10.
The prescribed final time is 0.02; 32 fixed auxiliary corrections separate
inverse error from spatial error. The new electron/chi field caps are separate
from retained original field caps. No existing accuracy gate is deferred for
this test.

Evidence directory:
[stage3-parallel-mms-4en7eogy](../build/validation/polarisation-drift-3d/stage3-parallel-mms-4en7eogy/).
It retains the generated input, criteria, runner and executable/library,
[hashes](../build/validation/polarisation-drift-3d/stage3-parallel-mms-4en7eogy/hashes.json),
[symbolic/input check](../build/validation/polarisation-drift-3d/stage3-parallel-mms-4en7eogy/symbolic-and-generated-files.log),
[Python check](../build/validation/polarisation-drift-3d/stage3-parallel-mms-4en7eogy/python-symbolic.log)
and [integrated log](../build/validation/polarisation-drift-3d/stage3-parallel-mms-4en7eogy/integrated.log).
Both serial reference levels completed at `t=0.02`. Before reusing them, the
resumed runner verified the executable hash, input hash and every grid/solver
override. The additional review also verifies that the retained libraries
match, the saved times are exactly `0, 0.005, 0.01, 0.015, 0.02`, and the
independently evaluated evolved-field errors agree with BOUT's native MMS
errors within the existing `2e-11` check.

The completed references are sufficient to evaluate the required 8-to-16
accuracy gates. They fail as follows:

| Assertion | Measured result | Required |
| --- | --- | --- |
| `order/16/Ni/linf` | `1.40129778792` | `>=1.5` |
| `order/16/Pe/linf` | `1.45645354090` | `>=1.5` |
| `order/16/Vort/l2` | `1.38409874406` | `>=1.7` |
| `order/16/Vort/linf` | `0.846937987918` | `>=1.5` |
| `order/16/Te/linf` | `1.44992919208` | `>=1.5` |
| `error/16/ion_chi` | L2 `0.00388090079533` | `<=0.001` |

All 18 recovery, continuity, auxiliary-residual, source-ownership and
nonzero-term checks pass for these two levels. The largest auxiliary residual
is `1.098594826220989e-11` against `1e-10`; source-ownership error is at most
`1.355252715606881e-20` against `1e-12`. All other field-error caps pass,
including vorticity L2 error `7.29546610942e-5` against
`8.165673925347495e-4`. Recovery tolerance and field-error caps are separate.

Checkpoint 3 is failed under the plan's first-failure policy. The serial
attempt had been interrupted after completing 16 cells to schedule remaining
work on eight Y ranks. That scheduling change predates accuracy evaluation;
the resumed 32-cell run is partial evidence only. Neither finer convergence,
the tighter 16-cell temporal repeat nor the alternate 8-cell decomposition
has completed. The 32-cell run was stopped with an interrupt to the identified
runner; its normal cleanup stopped the MPI launcher and retained the partial
output. See the [failure-interruption record](../build/validation/polarisation-drift-3d/stage3-parallel-mpi-ko75vl9f/checkpoint-failure-interruption.json).
No full three-level validation is claimed.

The vorticity order resembles the already recorded default-profile decline,
but that resemblance does not establish a cause or extend its deferral.
Chi L2 error decreases from `0.0138643205551` to `0.00388090079533`
(order `1.83691345622`) while exceeding its independent cap at 16 cells.
A small auxiliary residual establishes consistency with the numerical
equation, not accuracy against the continuum chi. Temporal separation and
fine-grid behaviour remain untested. The smallest proposed follow-up is to
complete those diagnostic comparisons, then localise the vorticity/chi errors
against the analytic reference, retaining all current accuracy thresholds.
No numerical-method change is made after this failed checkpoint.

The [checkpoint-3 decision](../build/validation/polarisation-drift-3d/stage3-parallel-mpi-ko75vl9f/checkpoint-3.json),
[completed-reference accuracy](../build/validation/polarisation-drift-3d/stage3-parallel-mpi-ko75vl9f/completed-reference-accuracy.json)
and [verified hashes](../build/validation/polarisation-drift-3d/stage3-parallel-mpi-ko75vl9f/verified-hashes.json)
retain all 66 evaluated assertions, norms and exact build/input identifiers.
The [resumed integrated log](../build/validation/polarisation-drift-3d/stage3-parallel-mpi-ko75vl9f/integrated.log)
records reference verification and the partial fine run. The original inputs,
outputs and scheduling-interruption record remain retained.

The [completed-reference evaluation script](../build/validation/polarisation-drift-3d/evaluate_stage3_references.py)
uses the frozen runner's physics checks and unchanged criteria; its
[log](../build/validation/polarisation-drift-3d/stage3-parallel-mpi-ko75vl9f/reference-evaluation.log)
records the accuracy decisions. Run it with:

```sh
.venv/bin/python build/validation/polarisation-drift-3d/evaluate_stage3_references.py
```

It returns status 1 for these six accuracy failures. The separate resumed-run
status `-2` records the deliberate interruption, not an evaluated fine-grid
accuracy result.

### Workflow correction and authorised electron-momentum diagnostics

The default [MMS runner](../tests/integrated/ion-energy-parallel-mms/runtest)
now evaluates each completed refinement pair immediately. It records all
checks in that block, then raises through `Study`'s result-writing context on
failure before launching another grid or dependent validation. Physics
failures also stop the next launch; a failed temporal comparison stops the
decomposition run. `Study.check()` retains its shared accumulating behaviour.
The runner has no accuracy-bypass option.

All six [workflow tests](../tests/python/test_parallel_mms_workflow.py) pass,
covering coarse physics, coarse/fine order, field-error and temporal failure,
plus successful sequencing. Two actual reuse runs save all 66 checks and the
same six numerical failures without launching any solver. See the
[workflow comparison](../build/validation/polarisation-drift-3d/stage3-electron-audit-ho56xmct/workflow-preservation.json),
[final reuse log](../build/validation/polarisation-drift-3d/stage3-electron-audit-ho56xmct/workflow-final-reuse.log)
and [test log](../build/validation/polarisation-drift-3d/stage3-electron-audit-ho56xmct/workflow-tests.log).
Formatting, Ruff and whitespace checks pass. Criteria are unchanged.

The maintainer separately authorised the electron-momentum RHS audit and
tighter 16-cell temporal comparison. These diagnostics keep checkpoint 3
failed and stages 4–5 pending. The costly 32-cell run has not been resumed.

The [read-only RHS audit](../build/validation/polarisation-drift-3d/audit_stage3_electron.py)
reproduces the saved-current decomposition at 16 cells. Current L2 error
grows from roundoff initially to `0.00512112803714`; the electron contribution
is `0.00512170900739`, while the combined ion contribution is
`2.88957090991e-5`. The reconstructed current-error vector agrees with the
saved current error to L2 `2.06711333940e-17`. Chi L2 error grows to
`0.00388090079533`. These are norms of vector contributions; their scalar
norms are not added to reconstruct the total.

The source audit finds the correct electron mass, charge and force signs:
`NVe=ne*ve/1836`, with electric force `+ne*Grad_par(phi)` and pressure force
`-Grad_par(Pe)`. The enabled continuum reference also contains parallel
inertia and ExB transport. BOUT adds the MMS source after the physical RHS;
the saved `ddt(NVe)` includes it. The generated source agrees with the
independent symbolic source to L2 `5.49275288882e-18` at the initial 16-cell
state. This comparison therefore subtracts the exact total derivative,
rather than interpreting the saved derivative as an unforced physical RHS.

At the initial state, the uniform identity-transform slab permits an exact
reconstruction of the numerical C2 parallel pressure and electric gradients.
Separating electric derivative truncation from the recovered-potential error
gives the following normalised force-density errors:

| Initial 16-cell electron RHS contribution | L2 error |
| --- | --- |
| Total saved RHS versus manufactured total derivative | `1.65095447491e-4` |
| Pressure-gradient truncation | `1.74997029524e-4` |
| Electric-gradient truncation on the analytic potential | `4.27551729198e-5` |
| Electric force from recovered-potential error | `8.22948582079e-5` |
| Vector sum of those three force errors | `1.65090766047e-4` |
| Remaining RHS error after subtracting their vector sum | `9.38541570388e-7` |

The remaining error is about `0.57%` of the total RHS-error L2 norm. The
initial total RHS error decreases from `6.57912755031e-4` at 8 cells to
`1.65095447491e-4` at 16 cells (order `1.99459594423`). This identifies the
initial force discretisation and potential recovery as the dominant measured
RHS discrepancy. It does not isolate which contribution controls the coupled
final-time current or vorticity error. Electron current multiplies momentum
error by `1836`, making these small force errors consequential. No operator
or manufactured forcing has been changed.

See the full [RHS audit JSON](../build/validation/polarisation-drift-3d/stage3-electron-audit-ho56xmct/electron-rhs-audit.json)
and [audit log](../build/validation/polarisation-drift-3d/stage3-electron-audit-ho56xmct/electron-rhs-audit.log).
The separate [temporal diagnostic](../build/validation/polarisation-drift-3d/diagnose_stage3_temporal.py)
uses the retained executable/library, identical input and serial decomposition,
with `atol=1e-14`, `rtol=1e-12`, `max_timestep=0.0005`. It checks the original
5%-of-spatial-error limits and absolute floor `1e-11`. Its initial sandbox
attempt failed MPI startup before any output. The first native attempt hit
the `600 s` no-saved-output watchdog at `t=0`, without reaching an evaluable
temporal comparison; that failed execution remains retained. A separate
[extended temporal diagnostic](../build/validation/polarisation-drift-3d/diagnose_stage3_temporal_extended.py)
used the same serial grid, input and numerical settings, an
`1800 s` output watchdog and `7200 s` total execution bound. These are
diagnostic execution bounds; default criteria and all accuracy limits are
unchanged.
The [diagnostic hashes](../build/validation/polarisation-drift-3d/stage3-electron-audit-ho56xmct/diagnostic-hashes.json)
and [native log](../build/validation/polarisation-drift-3d/stage3-electron-audit-ho56xmct/temporal-native.log)
retain the first diagnostic configuration. The
[extended log](../build/validation/polarisation-drift-3d/stage3-electron-audit-ho56xmct/temporal-native-extended.log)
and [extended hashes](../build/validation/polarisation-drift-3d/stage3-electron-audit-ho56xmct/temporal-extended/diagnostic-hashes.json)
retain the retry. The [audit source fingerprints](../build/validation/polarisation-drift-3d/stage3-electron-audit-ho56xmct/audit-source-hashes.json)
identify the generator, workflow code and inspected C++/BOUT sources. A diagnostic pass cannot advance
the failed spatial-accuracy checkpoint.

The extended attempt also hit its output watchdog. It ended after
`1805.557191584 s`, with the last saved output at `t=0`; its sole evaluated
check is the failed execution gate. The physics and 32 temporal-difference
checks were not reached. Both native attempts are preserved, including their
normal runner cleanup and failing exit status. No subsequent validation was
launched.

The [temporal diagnostic summary](../build/validation/polarisation-drift-3d/stage3-electron-audit-ho56xmct/temporal-summary.json)
and [extended raw result](../build/validation/polarisation-drift-3d/stage3-electron-audit-ho56xmct/temporal-extended/ion-energy-parallel-mms-m9wvphx2/result.json)
record this execution failure separately from the six spatial-accuracy
failures. The last saved time does not establish the internal integration
time reached before termination. These logs do not distinguish slow RHS
evaluation from small/rejected integration steps. The next diagnostic should
expose CVODE's internal step/iteration counts and RHS timings, with sufficient
output-progress visibility to complete the tighter 16-cell comparison before
prioritising 32 cells. The final-time current error's temporal contribution
remains unknown; the established initial force discrepancy is independent of
that missing comparison. No further numerical-method changes or runs follow
this failed bounded diagnostic without review.

## Constructor fix and review amendment

The conservative constructor formerly evaluated
`ion["laplacian"]["type"].as<std::string>()` eagerly as the argument to a drift
fallback. Before ion construction that option had no default value, including
when an explicit drift inverse type was supplied. The constructor now:

- Uses BOUT++'s `LaplaceFactory::default_type` if the ion type is omitted.
- Skips the ion-type read when the drift type is explicitly set.
- Preserves explicit ion solver types as the fallback.

`DriftBeforeIonResolvesDefaultAndExplicitInverse` directly exercises both
reported failures. `ReversedComponentConstructionPreservesRHS` constructs the
real scheduler with drift before ion energy/species, exercises a nonuniform
X-Z state with autosorting, and compares all evolved RHS values plus recovered
phi, chi, current and energy transport bitwise. The
[review probe](../build/validation/polarisation-drift-3d/stage1-review-efaf25jj/constructor-review-probe.log)
is preserved. The focused two-test run passed in 0.129 s and the final affected
80-test run passed in 1.527 s.

The [plan amendment](polarisation-drift-3d-analysis-plan.md#maintainer-amendment-after-stage-1-review)
explicitly defers `test_reynolds_wall_curvature_converges`, retains its
`>1.5` Reynolds threshold, and continues to prohibit Reynolds extraction or
changes to `vorticity.cxx`/`vorticity.hxx`. The previously measured boundary
orders `0.0252637567`, `0.0048517793` remain documented. New assertions still
follow the stop/report policy. A second maintainer amendment explicitly defers
only the eight existing default 8-to-16 MMS order gates named below, retaining
each `>=1.7` threshold and raw failure. Long-profile and new parallel-current
MMS gates are excluded.

## Stage-1 preservation evidence

The conservative extraction remains as described by the
[field contract](polarisation-drift-3d-field-contract.md): ion recovery owns
Reynolds/viscosity and publishes mandatory drift inputs; drift owns the
auxiliary inverse, fixed corrections, conservative face transport and uniform
wall products; finalisation consumes separately identified species sources.
Numerical diagnostic names/conversions are preserved. Version/owner metadata
changes are intentional; analysis now supports ion-energy versions 4, 5 and 6.

| Check | Acceptance | Current result |
| --- | --- | --- |
| Checkpoint 0 | Maintainer acceptance; kinetic exchange `>=5e-5` | Accepted; measured `8.533155513030926e-5` |
| Constructor order | Default/explicit drift solver works before ion construction | Passed |
| Affected C++/legacy tests | All pass | 80/80; vorticity and legacy fixtures unchanged |
| Retained C++ snapshots | Bitwise equality of valid numerical values | 67/67 snapshots agree |
| Closed/source default profiles | Existing criteria and bitwise preservation | Prior completed before/after runs pass; 121 field histories identical |
| Default MMS, original solver | Existing orders/errors/recovery thresholds | **Failed:** 49/57 checks pass; eight orders below 1.7 |
| MMS recovery/diagnostics | `1e-10` for recovery and diagnostic-state checks | All original recovery gates pass; diagnostic toggle exactly zero |
| MMS final-time L2 field errors | Separate normalised `criteria.json["errors"]` limits: Vort `8.165673925347495e-4`, Eion `1.1908840233405635e-3`; other fields retain their own limits | All original field-error gates pass |
| Completed MMS before/after cases | Bitwise physical-field and owned-face preservation | All five cases complete; 148 physical field histories and 234 numerical histories identical; all 57 assertion values identical |
| Wall original cases | Existing budgets, diagnostics, exchange and restart limits | All 13 executions complete; 204/204 checks pass |
| Wall migrated cases | Existing criteria and bitwise preservation | All 13 executions complete; 204/204 checks pass; 343 physical field histories and 566 numerical histories bitwise identical |
| MMS generator consistency | Symbolic checks and committed generated files agree | All four variants pass; regeneration changes section ordering only |
| Publication Python checks | All pass | Prior 85/85 pass; unchanged by constructor fix |
| Known Reynolds accuracy test | Retain test/threshold; report separately | Deferred by maintainer; existing failure remains enabled |
| Checkpoint 1 | Required preservation/default profiles pass, named accuracy exceptions only | **Passed under named existing-accuracy deferrals; raw failures and thresholds retained** |
| Checkpoint 2 | Parallel/sheath and preservation requirements | Passed with observed ordering roundoff accepted |
| Checkpoint 3 | New parallel-current MMS accuracy and sensitivity requirements | **Failed:** six new 8-to-16 accuracy gates; fine/temporal/MPI checks incomplete |
| Checkpoints 4–5 | Analysis/tokamak requirements | Not started following checkpoint-3 failure |

The optional unit recorder now checks whether `finally()` actually evaluated
`ddt(Eion)` before saving it. Allocation alone did not establish validity:
recovery-only tests could expose stale derivative memory. The comparison omits
only those undefined baseline entries; all evaluated derivatives remain
included, and no valid value differs. Earlier raw snapshots/logs are preserved.

MMS inputs were regenerated from the migrated base/generator after the
file-consistency check detected stale section ordering. A separate comparison
establishes identical settings and forcing expressions for every variant.
The runs use frozen before/after copies, so regeneration did not alter them.

## Existing default MMS accuracy failures

All reported orders use final-time L2 errors at 8 and 16 X/Z cells, with Y
resolution twice X/Z. Each run reaches `t=0.1` in 20 saved intervals of 0.005.
The original default profile took 1029.014 s on this validation environment.

| Variant / field | Historical recorded order | Fresh pre-extraction order | Required minimum |
| --- | ---: | ---: | ---: |
| baseline / Ni | 1.6650549299 | 1.6644897885 | 1.7 |
| baseline / Pe | 1.6731797065 | 1.6731816904 | 1.7 |
| baseline / Vort | 1.6903147855 | 1.3754612404 | 1.7 |
| baseline / Eion | 1.6201750290 | 1.6192537073 | 1.7 |
| curvature / Ni | 1.6943989984 | 1.6939794427 | 1.7 |
| curvature / Pe | 1.6144927754 | 1.6143863998 | 1.7 |
| curvature / Vort | 1.6912157236 | 1.3763096403 | 1.7 |
| curvature / Eion | 1.6667191125 | 1.6673157725 | 1.7 |
Other requested field orders pass. The retained solver/library and original
inputs predate the extraction; the resumed loader asserts dump version 4 for
every original case and version 5 for every migrated case. All five
paired MMS cases agree bitwise, including the fields exhibiting failed order
gates. This establishes preservation for the full default profile;
it does not establish the origin of every accuracy failure.

All eight historical orders in `criteria.json["observed"]` were already below
1.7. Six fresh results differ by less than 0.001, closely reproducing those
recorded limitations. Vorticity differs substantially: baseline falls from
1.6903147855 to 1.3754612404 and curvature from 1.6912157236 to 1.3763096403.
The older passing-default-suite narrative in
`docs/test-consolidation-validation.md` is inconsistent with those recorded
criteria. `observed_long` corroborates the historical 8-to-16 orders, but its
long-profile gates are not included in this deferral.

The criteria revision audit locates these historical values at commit
`fc5ecea0` (2026-11-09). That commit's validation notes explicitly record eight
order failures and a cadence comparison within `6.02e-13` relative L2; it did
not change the MMS forcing or thresholds. Since `734e7bfe`, the all-cell
`manufactured_errors` norm, continuum generator and base input are unchanged
apart from the stage-1 option-owner migration. Later pre-extraction ion-energy
commits changed the Reynolds construction, so a controlled historical-build
comparison is needed to locate the decline; this audit does not select a cause.
See the [criteria revision audit](../build/validation/polarisation-drift-3d/stage1-mms-policy-f7siip2t/historical-mms-criteria-revisions.json)
and [cadence commit](../build/validation/polarisation-drift-3d/stage1-mms-policy-f7siip2t/historical-mms-cadence-commit.txt).

The vorticity decline is a separate unresolved accuracy investigation. Both
fresh failures occur in the pre-extraction executable/library, and the completed
paired fields agree bitwise. This rules out the extraction as their trigger;
it does not establish which earlier code, dependency or measurement change
caused the decline. No forcing, boundary closure or accuracy limit was changed
to make them pass. The new parallel-current MMS must meet its own accuracy
requirements.

The retained pre-extraction vorticity errors are zero initially. Their final
maxima occur at physical X cells 3/8 and 7/16 (zero-based from the inner wall),
inside the domain. Excluding two physical cells at each radial wall gives L2
orders 1.6382686391 (baseline) and 1.6391324672 (curvature); final Linfinity
orders are 0.9287879655 and 0.9443873520. These diagnostic norms do not replace
the all-cell gates, and do not establish a Reynolds cause. See the
[radial errors](../build/validation/polarisation-drift-3d/stage1-mms-policy-f7siip2t/fresh-vorticity-radial-errors.json)
and [norm comparison](../build/validation/polarisation-drift-3d/stage1-mms-policy-f7siip2t/fresh-vorticity-norm-comparison.json).
An earlier progress update incorrectly described all baseline refinement
checks as passing after reading only the log tail. The full result JSON and
the table above supersede that statement; all eight failures remain enabled.

## Commands, build identity and retained artifacts

Current completion directory:
[stage1-mms-policy-f7siip2t](../build/validation/polarisation-drift-3d/stage1-mms-policy-f7siip2t/).
It preserves the report and policy before this amendment and contains the
[historical comparison](../build/validation/polarisation-drift-3d/stage1-mms-policy-f7siip2t/historical-mms-order-comparison.json).
Previous review directory:
[stage1-review-efaf25jj](../build/validation/polarisation-drift-3d/stage1-review-efaf25jj/).
Earlier evidence:
[stage1-vrgpz6qe](../build/validation/polarisation-drift-3d/stage1-vrgpz6qe/).
The earlier full report is retained as `report-before-constructor-fix.md`.

- Hermes: `779c4ff94f68e07fcfb725cc2fd576d9c1d76424` plus current changes.
- BOUT++: `97411814b32b23b2e8dbd1465e4160c072f9c5b3`; prior googletest change preserved.
- xhermes: `0a2f0325295a94f5a7567894917cb2bc9cef4c31`; no internal checkout changes.
- Python: `.venv/bin/python`, 3.13.3; installed editable xhermes.
- Build: `build/`, Release, CHECK=2, metric-3D OFF; ordinary 3D fields supported.
- All resumed integrated runs use one MPI rank, `FI_SOCKETS_IFACE=en0`, outside
  the sandbox because OFI initialisation is restricted inside it.
- Original runs set `DYLD_LIBRARY_PATH` to retained `original-build/` so the
  original executable resolves its original Hermes library. The dump-version
  assertion prevents accidental use of the rebuilt library.
- The required tokamak grid remains available, 36×48, SHA256
  `84c73a99d5b24f13d6215465c833902bf7ec229698fb0b4634bfd6afcec8d6cf`.
  No tokamak run was performed.

Build and focused tests:

```sh
cmake --build build --target hermes-3 hermes_unit_tests -j 4
FI_SOCKETS_IFACE=en0 build/hermes_unit_tests '--gtest_filter=*IonEnergy*:*PolarisationDrift*:*VorticityTest*'
.venv/bin/python tests/integrated/ion-energy-mms/mms.py --all --check --check-files
```

`run_preservation.py before|after mms|walls` invokes frozen copies of the actual
runners, with per-case timeouts 3600 s for MMS and 4800 s for walls. It reuses
only verified completed cases with matching executable/input hashes and
command-line overrides, and records their origins. The completion driver keeps
raw MMS failures and writes a separate `checkpoint-policy.json`; it does not
change the test result or exit status. All cases
retain `BOUT.inp`, resolved `BOUT.settings`, logs/dumps and `result.json`.
The interruption is separate from numerical failures in completed runs.
`run_wall_case.py` produces the two longest remaining independent migrated wall
cases. `run_preservation_parallel.py` waits for their completion, validates their
hashes/overrides and replays the full wall-profile audits. The interrupted serial
attempt is retained separately; it is not counted as a completed validation.

Key machine-readable evidence:

- [Checkpoint-1 decision](../build/validation/polarisation-drift-3d/stage1-mms-policy-f7siip2t/checkpoint-1.json).
- [Original wall result](../build/validation/polarisation-drift-3d/stage1-mms-policy-f7siip2t/before-runs/before-input-walls-vcpa5ho2/result.json),
  [migrated wall result](../build/validation/polarisation-drift-3d/stage1-mms-policy-f7siip2t/after-runs/after-input-walls-fb_q0sf7/result.json)
  and [all-assertion comparison](../build/validation/polarisation-drift-3d/stage1-mms-policy-f7siip2t/wall-checks-comparison.json).
- [Original MMS result](../build/validation/polarisation-drift-3d/stage1-review-efaf25jj/before-runs/before-input-mms-dby27zmo/result.json)
  and [failed assertions](../build/validation/polarisation-drift-3d/stage1-review-efaf25jj/mms-failed-checks.json).
- [Completed migrated MMS result](../build/validation/polarisation-drift-3d/stage1-mms-policy-f7siip2t/after-runs/after-input-mms-o0py0uld/result.json),
  [policy classification](../build/validation/polarisation-drift-3d/stage1-mms-policy-f7siip2t/after-runs/after-input-mms-o0py0uld/checkpoint-policy.json)
  and [all-assertion comparison](../build/validation/polarisation-drift-3d/stage1-mms-policy-f7siip2t/mms-checks-comparison.json).
- [Original MMS log](../build/validation/polarisation-drift-3d/stage1-review-efaf25jj/before-mms.log).
- [Current preservation comparison](../build/validation/polarisation-drift-3d/stage1-mms-policy-f7siip2t/preservation-comparison.json):
  completed cases, field counts, exact differences and incomplete case outcomes.
- [C++ log](../build/validation/polarisation-drift-3d/stage1-review-efaf25jj/unit-validated.log),
  [JUnit XML](../build/validation/polarisation-drift-3d/stage1-review-efaf25jj/unit-validated.xml)
  and [valid snapshot comparison](../build/validation/polarisation-drift-3d/stage1-review-efaf25jj/unit-validated-comparison.json).
- [Generated-input check](../build/validation/polarisation-drift-3d/stage1-review-efaf25jj/mms-generated-recheck.log)
  and [settings comparison](../build/validation/polarisation-drift-3d/stage1-review-efaf25jj/mms-canonical-settings-comparison.json).
- [Manifest](../build/validation/polarisation-drift-3d/stage1-review-efaf25jj/manifest.json),
  [build hashes](../build/validation/polarisation-drift-3d/stage1-review-efaf25jj/build-hashes.json)
  and [tracked patch](../build/validation/polarisation-drift-3d/stage1-review-efaf25jj/changes.patch).

The constructor source, tests and generated inputs pass formatting/whitespace
checks; Python lint passes. `vorticity` and the Reynolds discretisation remain
unchanged. No parallel-current/sheath equations were changed during the stage-1 runs;
stage-2 changes and their failed checkpoint are recorded above.

## Next step

Review the six new parallel-current MMS accuracy failures before resuming
dependent analysis/tokamak work. Stage-2 X-Z preservation has passed. Retain the
historical vorticity decline as a separate accuracy investigation with both
sets of measurements and their thresholds. The new parallel-current MMS must
pass its own accuracy requirements. Reynolds extraction and `vorticity`
changes remain outside this task.

The default-runner stop correction and initial electron RHS audit are
complete. Diagnose the tighter 16-cell execution limitation and complete its
temporal comparison before prioritising the costly 32-cell diagnostic.
