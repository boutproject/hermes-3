# Ion energy test follow-ups

Parallel-current divergence, `DivJextra` and `DivJcol` are now included in
the auxiliary and vorticity equations, with completed sheath current and
electrical-work diagnostics. The maintainer accepts the observed
`8.881784197001252e-16` trajectory-sorting difference as probably harmless
ordering roundoff; raw bitwise failures remain retained. Identical-state C++
RHS, target formulas, continuity, electrical work and serial/Y-MPI checks pass.
Stage-2 X-Z preservation also passes: 207 common numerical histories are
bitwise identical. Checkpoint 2 passes under this narrow acceptance. The new
independent-electron-momentum MMS is implemented but checkpoint 3 remains
failed on six new 8-to-16 accuracy gates. Its runner now checks each completed
refinement pair and writes failure results before launching more validation.
No existing accuracy deferral applies to it. The electron RHS audit identifies
parallel pressure/electric forces, including recovered-potential error, as the
dominant initial discrepancy. Tighter serial 16-cell diagnostics hit both the
600-second and extended 1800-second output watchdogs before saving a nonzero
time. Temporal sensitivity remains unassessed. Inspect internal CVODE step
counts and RHS timings before completing that comparison and prioritising
the costly 32-cell diagnostic; stages 4–5 remain pending.
The Y-varying `ion-energy-closed` input also exposed a `zero_current` closure
limitation: electron pressure work has no parallel force-balance/kinetic
counterpart. The revised strict X-Z benchmark uses Y-invariant fields, and
analysis reports unpaired work separately from heating. Its energy-conservation
check passes; the maintainer accepted checkpoint 0 with the parallel kinetic
exchange minimum revised to `5e-5` of initial energy. See
[the polarisation-drift review](polarisation-drift-3d-analysis-report.md).
Conservative drift has been extracted with bitwise-preserved completed X-Z
outputs, and the constructor-order defect is fixed. The maintainer deferred
the known Reynolds wall-accuracy failure while retaining its test and `1.5`
threshold. Checkpoint 1 subsequently found eight original default MMS order
failures against `1.7`. The maintainer now defers only those eight named
default 8-to-16 accuracy gates, retaining tests and thresholds; checkpoint 1 preservation is complete. All five MMS and 13 wall cases
agree bitwise with the originals; both wall profiles pass all 204 assertions. All eight historical orders were already
below 1.7; six fresh values reproduce them closely. Baseline/curvature vorticity
fall from about 1.69 to 1.38 before extraction and require a separate accuracy
investigation. The old passing-suite narrative conflicts with the recorded
criteria. Long-profile and new parallel-current MMS gates remain outside this
deferral. The parallel-current/sheath implementation, initial failed checkpoint
and subsequent maintainer acceptance are documented in the review report. For the separate vorticity decline, compare frozen MMS inputs and
error norms with historical ion-energy builds at a fixed BOUT++ revision;
locate the change before proposing a numerical fix. Fresh error maxima lie
inside the radial domain, and excluding two wall cells still gives orders
about 1.638/1.639. This evidence does not attribute the decline to the known
Reynolds wall failure. Reynolds extraction and `vorticity` changes remain separate.
The removed cases are preserved at commit `3831b4d`:

- [Open parallel transport](https://github.com/boutproject/hermes-3/tree/3831b4d/tests/integrated/2D-ion-energy-boundary):
  prescribed ion flow through open Y boundaries.
- [Electron parallel energy and pressure work](https://github.com/boutproject/hermes-3/tree/3831b4d/tests/integrated/2D-ion-energy-parallel):
  restore with current closure. Its numerical heating audit found
  `FV::Div_par_fvv_heating` using `n_edge*(a +/- v_edge/2)` where the split
  kinetic balance needs `n_edge*(a +/- v_edge)/2`. Check this against the
  `operator-fixes` BOUT++ branch before restoring the gate. The current
  checkout (`97411814`) already uses the corrected split particle flux; the
  analysis library and its kinetic split identity now use that form.
- [Sheath and tokamak relaxation](https://github.com/boutproject/hermes-3/tree/3831b4d/tests/integrated/tokamak-ion-energy-relax-boundary):
  COMPASS geometry, closed-core surface averaging, two targets and six ranks.
  Restore with the `long` label.
- Add a nonzero-parallel-current variant to `ion-energy-mms`. Varying magnetic
  field and the alternative Reynolds stress split remain separate work.

Reusable functions in the removed
[ion_energy_tokamak_test.py](https://github.com/boutproject/hermes-3/blob/3831b4d/tests/integrated/ion_energy_tokamak_test.py)
are `read_ranks` (topology-aware rank data), `boundary_rates` (wall histories),
`ode_error` (boundary ODE integration), and `sheath_state` (target reconstruction).
Port these into `energy_analysis` when the corresponding physics is supported.

MMS source ownership guards in `evolve_energy`, `neutral_mixed` and
`neutral_full_velocity` are deferred. Density and pressure guards are retained
and tested; other MMS variable sources need the same ownership rule before
those components are used in expression-source MMS cases.

The consolidated source test retains collisional heat exchange but omits
`braginskii_friction`: with `zero_current`, electrons have no evolved kinetic
reservoir to receive their frictional mechanical work. Heat exchange cancels
between species; friction does not conserve the chosen plasma reservoir under
that closure. `ParticleInjectionConservesTransferredEnergy` and
`ParticleRemovalCarriesItsKineticEnergy` cover momentum-source conversion.
