# Relaxing potential boundaries in the total ion energy model

Status: implemented. The production component uses an auxiliary polarisation
potential `chi` and solver-managed phi midpoint state stored in physical
radial boundary cells of `Vort`. The sections below record the design contract.

Validation includes C++ storage/packing, recovery, restart and mixed-metric
face tests; independent Python reconstruction tests; the separate slab
diagnostic/refinement/MPI/restart matrix; and a six-rank COMPASS/sheath matrix.
The COMPASS audit checks independent sheath source assembly and parallel ion
power, including core/private-flux interface cancellation. It explicitly
retains incomplete coverage for general tokamak electrical work and electron
exchange rather than claiming a complete independent tokamak power balance.

Enable `phi_boundary_relax` in
[`evolve_vorticity_ion_energy.cxx`](../../src/evolve_vorticity_ion_energy.cxx),
using the relaxation law and core averaging in
[`vorticity.cxx`](../../src/vorticity.cxx). Retain polarisation current in the
vorticity equation. Close the normal polarisation fluxes of particles, parallel
momentum, enthalpy, parallel kinetic energy and electrical work at physical
boundaries. These closures apply to faces, including the `g13` contribution;
they do not zero the RHS of boundary-adjacent cells.

## Starting point before implementation

- The constructor reads `phi_boundary_relax`, `phi_boundary_timescale` and
  `phi_core_averagey`, and sets `INVERT_SET` radial solver flags. The actual
  relaxation RHS is missing. The proposed solver-managed approach does not
  need simulation time or a last-update clock.
- Recovery iterates on `Phi = phi + psi`, with separate pressure boundary
  conditions and the finite-volume forward operator `Div_a_Grad_perp_XZ`.
  Homogeneous correction boundaries alone do not preserve a prescribed `phi`
  boundary when the pressure potential `psi` changes between iterations.
- The polarisation potential solve uses the full forward operator, then
  unconditionally overwrites its guards with `neumann_o2`. This can invalidate
  the residual checked immediately before the overwrite.
- `divPol(weight)` already passes `bndry_flux=false` to every constituent
  transport operator. `Div_a_Grad_perp_XZ` closes the entire radial face flux,
  including `g13`; XPPM closes physical X and, when enabled, open Y faces.
  Particle and momentum sources and the combined `DivEpol` use this helper.
- `DivJpol` currently diagnoses this closed transport divergence. The Python
  flux audit rejects potential relaxation and reproduces the old post-solve
  Neumann rule.

## 1. Evolve the phi boundary using Vort boundary cells

Use BOUT's `Vort:evolve_bndry=true` support to carry the boundary midpoint
potential `b_phi(y)` as solver state. This is preferable to advancing history
inside `transform`: boundary relaxation then participates in the same
adaptive stages, rejected-step rollback, error control and restart mechanism
as the interior variables. For a fixed target the ODE gives the same
exponential relaxation law; a changing target is evaluated from the recovered
potential at the current solver stage.

The relevant implementation is in
[`Solver::add`, `pre_rhs`, `post_rhs` and variable packing](../../external/BOUT-dev/src/solver/solver.cxx).
The 3D registration reads `evolve_bndry` before storing the variable metadata.
Set the option before `solver->add(Vort, "Vort")`; the constructor currently
registers `Vort` before reading the relaxation option and must be reordered.
Enable this specifically for `Vort`, rather than through `all:evolve_bndry`.
Validate a finite, positive normalised relaxation timescale and the required
radial guard width.

### Solver storage and physical vorticity

| Part of the solver's Vort field | Meaning in relaxation mode |
| --- | --- |
| `RGN_NOBNDRY` | Physical vorticity, evolved with the existing vorticity equation. |
| Physical `RGN_INNER_X` / `RGN_OUTER_X`, at physical Y indices | Boundary midpoint potential in normalised potential units. |
| Other `RGN_BNDRY` cells, including unused open-Y/corner cells | Explicitly initialised unused storage, with zero RHS. |
| Processor and periodic halos | Ordinary vorticity halo data supplied by communication. |

`evolve_bndry` includes the full physical `RGN_BNDRY`, not just the first
radial layer or just X boundaries. Account for every included cell. A simple
representation replicates `b_phi(y)` across Z and radial boundary layers,
initialises those copies identically, and evolves each copy with the same
target and its own stored value. Decode the Z mean in the first radial layer
to enforce a Z-independent boundary. Identical copies should remain identical;
test this invariant. Do not project or overwrite evolved storage in
`transform`. Communication must preserve physical boundary storage while
supplying processor/periodic halos.

Keep a separate physical-vorticity field for calculations and publication to
`state["fields"]["vorticity"]`. Copy the physical interior and reconstruct its
ordinary vorticity boundary/halo values using separate physical boundary rules.
In particular, XPPM reads neighbouring vorticity values when limiting the flux
at an interior face next to a boundary, even with boundary transport disabled.
Using the stored potential values there would mix two different physical
quantities and change interior Reynolds/polarisation transport. Use this
physical field consistently in the Reynolds and polarisation helpers.

### Boundary operators and recovery

BOUT applies value boundary operators in `pre_rhs` and derivative operators
in `post_rhs` even when `evolve_bndry=true`. Install no-op boundary operators
(`none`/`free`) on the raw storage field so these calls preserve both the
stored potentials and their relaxation rates. Preserve the physical-vorticity
boundary configuration separately; copying the raw field's no-op boundary
metadata is insufficient. Check effective per-region settings as well as
`bndry_all`, including corners and open-Y boundaries.

In `transform`, decode `b_phi` from the solver state and hold it fixed during
pressure/potential recovery. Reimpose

\[
\phi_g=2b_\phi-\phi_i
\]

before constructing `Phi`, after every potential correction, and after each
pressure/`psi` update. The stored quantity is the physical boundary midpoint,
not the reconstructed guard value `phi_g`. Populate the second phi guard
consistently with the legacy rule and fill Y halos without overwriting radial
conditions. Use homogeneous correction solves and enforce the nonhomogeneous
boundary on the total iterate. Assess the full forward residual and prescribed
phi midpoint after the final boundary/communication operations.

In `finally`, initialise all of `ddt(Vort)` to zero, assign the existing
vorticity RHS only in physical interior cells, then assign each stored radial
boundary rate as

\[
\dot b_\phi=\frac{a(\phi)-b_\phi}{\tau}.
\]

The target `a` is the adjacent interior Z average, with an optional
closed-core Y average at the inner boundary. Use the converged potential of
this stage. Reduce both sum and actual physical sample count for core
averaging; the legacy division by `LocalNy` includes Y guards. The outer
boundary target is independent for each physical Y index. This RHS does not
advance or mutate the boundary state during evaluation. Remove the last-update
clock and history-update logic from the planned implementation.

### Initialisation, restart, output and solver tolerances

On a fresh start, initialise the stored midpoint from the initial phi boundary
before the solver first packs its state. The vorticity initial-profile
expression in guard cells is not a potential initial condition. New restart
files save the encoded `Vort` field through the solver, including its boundary
state. Add a small static encoding/version marker and validate it when loading.
For older restarts, initialise boundary storage once from restored `phi`;
repeated restart/output hooks must never overwrite an evolved boundary.

Keep ordinary output `Vort` as physical vorticity with physical guard values,
and keep the encoded solver field in restart output. BOUT has separate normal
and restart output hooks, so those representations can be kept distinct.
Existing saved phi guards expose the boundary midpoint for independent
analysis without an extra volume output. Document the mixed-unit restart
encoding explicitly; consumers must not apply vorticity unit conversion to
its potential boundary cells.

The existing per-variable `Vort` absolute tolerance now controls both
normalised vorticity and normalised boundary potential. Verify that it is
adequate for both and include relaxation-timescale/time-solver refinement.
The boundary equation also adds a `-1/tau` contribution to the Jacobian;
include it in any applicable user preconditioner and test solver packing,
boundary operators, adaptive stepping and uninterrupted-versus-restarted runs.

## 2. Define the polarisation potential boundary condition explicitly

Let `L(f) = Div_a_Grad_perp_XZ(c, f)` with full boundary fluxes and fixed
`c = average_atomic_mass/B^2`. Retain

\[
L(\Phi)=\mathrm{Vort},\qquad
\dot{\mathrm{Vort}}=\mathrm{DivJReynolds}+\mathrm{DivJdia}.
\]

Do not close `L(Phi)` or add a second source to represent current from the
relaxing boundary. The full operator already includes that current. Retain the
existing boundary choices of each Reynolds-current term; opening all of their
faces would be an additional model change.

Use the agreed auxiliary-current contract for the field currently named
`dPhidt`: solve `L(chi) = dot(Vort)` in physical interior cells with
explicit homogeneous midpoint Dirichlet conditions for `chi` in relaxation
mode. This provides a nonsingular full-current solve using the existing
preconditioner and avoids adding a coupled RHS assembly stage. It is an
explicit approximation to the spatial polarisation current, including its
interior transport; it is not an exact derivative of the recovered potential.
These fluxes are formally small, and the approximation does not affect
conservation of the evolved total energy with the specified face closures.
The `dot(Vort)` used here is the physical vorticity RHS, excluding the
potential relaxation rates stored in its boundary cells.

Rename/document the internal field accordingly and give its diagnostic an unambiguous name, such as
`ion_polarisation_chi`, retaining potential-per-time units. If the existing
`ddt(Phi)` name is retained for compatibility, its metadata and analysis must
identify the approximation.

In this path, impose the `chi` boundary during defect correction and remove the
post-solve Neumann overwrite in relaxation mode. Check `L(chi) = dot(Vort)`
after all boundary operations. Close transport using face masks below; a
Neumann guard rule is unnecessary for that purpose.

The relaxing boundary rate therefore does not prescribe the auxiliary `chi`
boundary. An exact `dot(Phi)` would instead require `dot(b_phi) + dot(psi_b)`
and coupled density/momentum/energy rates; that derivative problem is outside
this implementation.

## 3. Keep current continuity and closed transport distinct

| Quantity | Physical boundary policy |
| --- | --- |
| Forward vorticity operator and polarisation current-continuity solve | Keep their full normal current fluxes. |
| Species particle flux `(A_s n_s / (Z_s rho)) * Jpol` | Zero normal polarisation transport faces. |
| Species parallel momentum flux `(A_s M_s / (Z_s rho)) * Jpol` | Apply the same closure. |
| Enthalpy flux `sum_s[(A_s/(Z_s rho)) * (5/2) p_s] * Jpol` | Apply the same closure. |
| Parallel kinetic energy flux `sum_s[(A_s/(Z_s rho)) * Kparallel_s] * Jpol` | Apply the same closure. |
| Polarisation electrical work `phi * Jpol` | Apply the same closure. |

Retain the common weighted face transport in `divPol` and its
`E_adv = phi + sum_s[A_s/(Z_s rho) * ((5/2)*p_s + Kparallel_s)]` coefficient.
Make the closure explicit in helper naming/comments and use one physical-face
policy for every constituent term. Test each energy contribution separately
so cancellation in the combined `DivEpol` cannot conceal an open flux.
Provide test access to the shared face fluxes through a small helper rather
than checking only a volume integral, where opposite boundary contributions
could cancel.

Retain `polarisation_sources` bookkeeping: their local parallel kinetic-energy
conversion is removed from the external source sum because their energy
transport is already included in `DivEpol`. Preserve tangential fluxes,
processor-interface transport and periodic connections. Preserve imposed
sources and sheath/parallel fluxes. In particular, `DivEelectric` also contains
diamagnetic and parallel electrical work and must not be zeroed wholesale.

Keep `DivJpol`'s existing meaning as the closed transport divergence, with clear
metadata. If needed, add optional full-current face diagnostics with separate
names. Species charge transport must sum to the closed transport current above
the density floor. At a physical boundary the retained current and transported
species charge intentionally differ under the requested model prescription.
Do not assert equality there or zero the whole boundary-adjacent divergence.

Audit remaining energy face terms, especially `DivEReynolds`, independently.
Closing the five requested polarisation contributions does not imply that
every existing energy boundary flux vanishes. Report any remaining flux and
its physical origin; extending closure to other flux families is a separate
boundary-model decision. Add no compensating boundary heating source.

## 4. Tests and energy analysis

Extend [`test_evolve_vorticity_ion_energy.cxx`](../unit/test_evolve_vorticity_ion_energy.cxx)
with the following checks:

- Fixed-target exponential relaxation on both radial boundaries, normalised
  rates, repeated/backtracking RHS calls that preserve solver state, adaptive
  integration and restart. Check equality of replicated boundary states.
- Real solver packing/unpacking and pre-/post-RHS boundary operations preserve
  potential state and rates. All evolved boundary cells have finite,
  explicitly assigned derivatives, including unused open-Y/corner cells.
- Changing stored boundary potentials cannot contaminate physical vorticity
  guards used by flux limiters. Compare reconstructed physical stencils with
  a separately manufactured vorticity field.
- Z averages and closed-core Y averages using physical sample counts, including
  singleton Y and an MPI check of core averaging.
- Recovery with changing pressure and unequal ion mass/charge: prescribed
  `phi` midpoint and full forward residual survive all iterations and halos.
- Manufactured full-current solve with a nonzero normal boundary current and
  nonzero `g13`, using the selected auxiliary boundary condition and excluding
  encoded potential rates from the physical current RHS.
- Nonzero unmasked current but exactly zero normal particle, momentum,
  enthalpy, parallel kinetic and polarisation electrical face fluxes. Check
  shared species charge partition and discrete Gauss balance, including
  boundary-adjacent cells and processor/tangential faces.
- The `phi_boundary_relax=false` regression suite and diagnostic on/off
  equivalence.

Add a separate `2D-ion-energy-relax-boundary/` integrated directory with its own
README, input, criteria and `runtest`. Use a nonperiodic radial slab, nontrivial
Z structure, initial boundary/interior mismatch and two ions with unequal
mass-to-charge ratios. Require measurable relaxation/current and independent
face-flux/source budgets. Run diagnostics on/off, restart and time-step/output
cadence refinement; add a radial MPI split. The existing
`2D-ion-energy-boundary` input has periodic X and imposed parallel flow and does
not exercise this feature.

Extend `tools/energy_analysis/fluxes.py` with the selected `chi` solve boundary
contract before removing its relaxation rejection. Reconstruct `psi` and
midpoint `phi` from existing fields/guards and use `bd.options`. For the
recommended auxiliary solve, the known homogeneous `chi` boundary allows an
independent sparse reconstruction without boundary-history outputs. Regular
saved-output differences of `Phi` remain cadence-dependent estimates of the
actual derivative and must not be compared as exact auxiliary `chi` values.

The boundary relaxation rate is reconstructible from the saved midpoint,
adjacent potential and configured timescale/averaging. No legacy boundary
history is needed. Keep detailed current maps optional. Record all remaining
source and boundary coverage explicitly, including Reynolds/electrical powers.

After the slab contracts pass, add a separate COMPASS/sheath case using six MPI
ranks as described in the [energy analysis plan](ion-energy-analysis-plan.md).
Check radial relaxation, zero requested polarisation transport and unchanged
source/sheath flux assembly together, including physical target faces and
internal topology connections. Any unsupported flux retains incomplete audit
coverage rather than an assumed zero contribution.

## Implementation order and completion criteria

1. Register solver-managed Vort boundary storage with no-op boundary operators
   and a separate physical-vorticity field.
2. Add boundary initialisation, RHS, restart handling and recovery enforcement.
3. Implement the chosen full-current solve boundary and explicit transport
   policy; clarify diagnostics.
4. Add unit/MPI tests, Python reconstruction and the separate slab test.
5. Extend to the tokamak/sheath boundary test with independently audited fluxes.

Completion requires converged full-current residuals, correctly relaxing
`phi` midpoints, zero normal flux for each requested polarisation transport
family, correct current/species distinction, reproducible time/restart
behaviour for solver-managed relaxation, and budgets with explicit source/boundary
coverage. A small global energy drift alone is insufficient evidence.
