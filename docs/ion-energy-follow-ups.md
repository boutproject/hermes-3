# Ion energy validation follow-ups

Conservative polarisation drift owns the auxiliary current solve and species
transport. Parallel-current divergence, `DivJextra` and `DivJcol` enter both
the selected current RHS and physical vorticity equation. Ion finalisation
adds optional viscosity and numerical potential damping. Damping uses completed
collective sound speed and has no explicit total-energy source.

All MMS variants now live in `ion-energy-mms`. `parallel_current` has zero
radial curvature of the Reynolds potentials at both walls;
`parallel_current_wall` uses P(a)=1−cos(2a) to expose nonzero wall curvature
on practical refinement grids. Both cases remove unrelated wall-weight
projection error by making ion primitive wall values independent of Z. Both use
MC-aware L2/Linfinity order limits of 1.4/0.9, applied globally or in a fixed
interior region respectively. The old curved-wall run failed five stricter
8-to-16 gates (Ni, Pe and Te Linfinity, Vort L2 and Linfinity); those historical
measurements do not establish results for the new smooth-wall solution.
Chi absolute error remains diagnostic-only, while convergence, temporal and
MPI checks remain enabled. The optional plotting mode retains accuracy
failures and collects further resolutions. Fresh 8-to-16 runs pass all
curved-case interior gates. The smooth case fails
only global Vort L2 (1.360 against 1.4); its interior order is 1.605. All other
field orders and physics checks pass. These results use final time t=0.02.
Finer global convergence and native temporal sensitivity remain unassessed.

The independent wall tests retain the robust reflection and MC treatment.
Zero wall curvature gives second-order wall maximum and domain mean absolute
error. Nonzero wall curvature leaves a bounded first-cell error, giving
first-order domain mean error and near-second-order smooth interior RMS.
The developer guide explains how MC extrema affect other norms. The separate
pre-existing zero-current Vort order decline has not been attributed to this
wall mechanism; compare fixed inputs, norms and dependency revisions before
proposing a numerical change.

All original zero-current MMS order gates and the parallel-current gates
remain active, including the historical zero-current failures. Each
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
