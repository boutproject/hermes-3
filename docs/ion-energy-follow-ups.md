# Ion energy validation follow-ups

Conservative polarisation drift owns the auxiliary current solve and species
transport. Parallel-current divergence, `DivJextra` and `DivJcol` enter both
the selected current RHS and physical vorticity equation. Ion finalisation
adds optional viscosity and numerical potential damping. Damping uses completed
collective sound speed and has no explicit total-energy source.

The separate `ion-energy-parallel-mms` currently fails five 8-to-16 order
gates: Ni, Pe and Te Linfinity, plus Vort L2 and Linfinity. Chi meets its spatial
orders; absolute chi error is diagnostic-only, with reference and relative RMS
saved. The runner writes results and stops before further cases after a failed
gate. Coupled temporal sensitivity remains unassessed. Initial production
electron-force probes show pressure/electric derivative and recovered-potential
errors, while small momentum-only temporal refinements are negligible compared
with their spatial errors. Those probes do not establish a final-time cause.

The independent wall tests retain the robust reflection and MC treatment.
Zero wall curvature gives second-order wall maximum and domain mean absolute
error. Nonzero wall curvature leaves a bounded first-cell error, giving
first-order domain mean error and near-second-order smooth interior RMS.
The developer guide explains how MC extrema affect other norms. The separate
pre-existing zero-current Vort order decline has not been attributed to this
wall mechanism; compare fixed inputs, norms and dependency revisions before
proposing a numerical change.

Existing default MMS acceptance excludes only its eight named 8-to-16 order
failures. Long-profile and parallel-current MMS gates remain active. Each
runner retains raw errors, orders and effective limits in `result.json`.
Conservation, source ownership and current identities do not establish passing
coupled accuracy.

The Y-varying `zero_current` closure has no independent electron kinetic
reservoir or parallel force balance. A strict closed-energy case must eliminate
unpaired electron pressure/electric work; the closed X–Z benchmark uses
Y-invariant fields. Analysis exposes this work separately from numerical
heating. Target formulas, continuity, electrical work and serial/Y-MPI
comparisons remain tested. Near-roundoff component-reordering trajectory
comparisons were removed.

The removed cases are available at commit `3831b4d`:

- [Open parallel transport](https://github.com/boutproject/hermes-3/tree/3831b4d/tests/integrated/2D-ion-energy-boundary):
  prescribed ion flow through open Y boundaries.
- [Electron parallel energy and pressure work](https://github.com/boutproject/hermes-3/tree/3831b4d/tests/integrated/2D-ion-energy-parallel):
  restore with current closure. Verify the corrected split particle flux in
  `FV::Div_par_fvv_heating` against the kinetic balance before restoring a gate.
- [Sheath and tokamak relaxation](https://github.com/boutproject/hermes-3/tree/3831b4d/tests/integrated/tokamak-ion-energy-relax-boundary):
  COMPASS geometry, closed-core averaging, two targets and six ranks.

The removed
[ion_energy_tokamak_test.py](https://github.com/boutproject/hermes-3/blob/3831b4d/tests/integrated/ion_energy_tokamak_test.py)
has topology-aware loading, boundary histories, wall-ODE and sheath helpers.
Use the existing energy-analysis adapters when extending geometry support.

MMS source ownership in `evolve_energy`, `neutral_mixed` and
`neutral_full_velocity` still needs the same injection checks used for
density/pressure before expression-source MMS can cover those components.
The consolidated source test omits friction under `zero_current`, because
there is no electron kinetic counterpart to its mechanical work. Collisional
heat exchange and moving-ion source conversion remain tested.
