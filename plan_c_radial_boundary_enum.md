# Plan C: one enum for physically consistent radial boundary conditions

## 1. Context

### 1.1 Component and prerequisites

- Component: `evolve_vorticity_ion_energy` in Hermes-3 (BOUT++ `next`, C++20).
  - `include/evolve_vorticity_ion_energy.hxx`
  - `src/evolve_vorticity_ion_energy.cxx`
- The component is new and not used in production. **No backward
  compatibility, deprecation warnings or option aliases.** BOUT++ rejects
  unused input options, so read each option only where it applies.
- **Prerequisite: Plan A** (`plan_a_phi_relaxation_velocity.md`) is
  implemented. This plan uses three things from it:
  - `phi_boundary_velocity` (κ), and the second-order wall law
    dφ_b/dt = κ̂ |∇x|_b D_in, computed by `computeWallRates()` in
    `transform_impl()` and written to `ddt(Vort)` by `writeWallRates()` in
    `finally()`.
  - The flux-surface averaging of the stored wall value and its rate at
    closed flux surfaces.
  - `applyPolarisationBoundary(chi, dc_inner, dc_outer)`, which takes DC wall
    values for χ.

  If Plan A is not yet implemented, implement its §3 first. It is
  independent of this plan.
- Keep the warm-start behaviour: the RHS must stay a deterministic function
  of the state. `phi = copy(phi_guess)` at the start of `transform_impl` is
  followed by `applyPotentialBoundary()`, and `acceptWarmStart()` runs only
  between steps.

### 1.2 Current radial boundary controls

All of these must be replaced or derived from the new enums.

| Item | Where | Current behaviour |
|---|---|---|
| `phi_boundary_relax` | constructor | Switches between two different φ boundary mechanisms |
| `phi_boundary_velocity` (after Plan A) | constructor | Relaxation velocity |
| `phi_core_averagey` | constructor | Flux-surface treatment at the core |
| `neumann_boundary_average_z` | `applyPressureBoundary()` | If true: p̃ = 0 at the midpoint, Neumann for the z-mean. If false: plain Neumann |
| `Vort:evolve_bndry` | constructor | Forced true only with relaxation |
| `Vort:bndry_*` (user input) | swap into `vorticity` | X and Y boundary conditions of the physical vorticity used in the XPPM stencils |
| `"neumann"` hard-coded | `calculateViscosity()` | Vorticity boundary for viscosity, separate from the above |
| `laplacian:inner/outer_boundary_flags` | `phiSolver` | Forced `INVERT_SET` with relaxation; otherwise user or Laplacian default (implicit Dirichlet, preserving the initial midpoint through homogeneous corrections) |
| `polarisation_laplacian:inner/outer_boundary_flags` | `chiSolver` | Validated to be `INVERT_AC_GRAD` (DC midpoint Dirichlet, AC Neumann) |
| `vEdotGradPsi.applyBoundary("neumann_o2")`, `DelpPhi_2B2.applyBoundary("neumann_o2")` | Reynolds section | Boundary conditions on derived quantities |
| local `Field3D chi = phi + psi; chi.applyBoundary("neumann_o2")` | `DivEReynolds` block | Shadows member `chi` |
| `Div_a_Grad_perp_XZ(...)` default `bndry_flux = true` | Reynolds Q term in `DivJReynolds`, `DivEReynolds`, viscosity | Wall fluxes partly open (g13 corner terms leak), while the XPPM calls pass `false` |
| `bndry_flux = false` | `divPol`, `DivEExB`, Reynolds XPPM | Closed |

The header documentation also lists options that are not read: `bndry_flux`,
`evolve_log`, `hyper_z`, `precon`, `E<name>:source` and
`source_only_in_core`. It also gives the `poloidal_flows` default as true;
the code default is false.

### 1.3 Physics of an impermeable radial wall

These conditions are common to every mode.

- **Potential.** φ uniform in z along the wall (φ̃ = 0 at the midpoint), so
  v_E·n = ∂_zφ = 0.
- **Pressure potential.** ψ uniform in z along the wall (ψ̃ = 0 at the
  midpoint), so U·n = ∂_zΦ = 0, where U = v_E + v_dia.
- **Consequences:**
  - No E×B or diamagnetic advective flux crosses the wall.
  - No momentum flux crosses the wall: the Reynolds and diamagnetic stresses
    carry a factor v_E^x or U^x.
  - No K⊥ energy flux crosses the wall.
  - Q = v_E·∇Φ = 0 on the wall.
- **Vorticity flux is not automatically zero.** In the split form used here,
  the natural wall flux is ½ c g^{xx}(∂ₓ∂_zφ ∂ₓψ − ∂ₓφ ∂ₓ∂_zψ). It vanishes
  for ψ = 0. Closing the Reynolds wall flux drops only this diamagnetic part,
  which is a statement that the wall is insulating for it.

The modes differ only in (1) how the mean wall potential evolves and (2) the
wall vorticity condition.

## 2. Goal and acceptance criteria

1. **New options.** The only radial boundary options are:
   - `radial_boundary`: `fixed` (default), `free_slip` or `no_slip`.
   - `phi_boundary_velocity` [m/s]: read only, and required, when
     `radial_boundary = no_slip`.
   - `polarisation_ac_boundary`: `neumann` (default) or `dirichlet`.

   Both enums are read with `BOUT_ENUM_CLASS_NS`.
2. **Removed options.** `phi_boundary_relax`, `phi_core_averagey` and
   `neumann_boundary_average_z` are removed.
3. **Rejected user settings.** If any of these are set in the input, throw a
   `BoutException` that names `radial_boundary`:
   - `laplacian:inner_boundary_flags`, `laplacian:outer_boundary_flags`
   - `polarisation_laplacian:inner_boundary_flags`,
     `polarisation_laplacian:outer_boundary_flags`
   - `Vort:bndry_xin`, `Vort:bndry_xout`
   - `Vort:evolve_bndry = false`
4. **One code path.** Every run with non-periodic X uses the stored-wall-state
   mechanism: `Vort` radial guard cells hold the φ midpoint, and
   `evolve_bndry` is forced. The old implicit-boundary path is deleted.
5. **Single source of truth.** One struct, filled in the constructor, holds
   every derived boundary setting. Every call site reads from it. A summary is
   printed to `output_info` at startup.
6. **Exact wall budgets.** Every wall flux, for charge, density, momentum,
   energy and K⊥, is zero by construction except those listed in §3.3. Tests
   in §5 check this from diagnostics.
7. **Clean build.** It compiles with `-Wall -Wextra` without new warnings, and
   the existing unit tests pass.

## 3. Design

### 3.1 Enums (header, namespace scope)

```cpp
#include <bout/bout_enum_class.hxx>
BOUT_ENUM_CLASS_NS(bout::ion_energy, RadialBoundary, fixed, free_slip, no_slip);
BOUT_ENUM_CLASS_NS(bout::ion_energy, PolarisationACBoundary, neumann, dirichlet);
```

Read them with:

```cpp
radial_boundary = options["radial_boundary"]
    .doc("Radial wall model: fixed, free_slip or no_slip")
    .withDefault(bout::ion_energy::RadialBoundary::fixed);
```

This pattern was checked to compile against BOUT++ `next`. Input strings are
lowercase.

### 3.2 Derived settings struct

```cpp
struct RadialBoundaryModel {
  bool walls{false};               // !mesh->periodicX
  RadialBoundary mode;             // from options
  PolarisationACBoundary chi_ac;   // from options
  bool evolve_wall_potential;      // mode == no_slip
  bool vorticity_wall_zero;        // mode == free_slip  (ω = 0 at midpoint)
  bool viscous_wall_flux;          // mode == free_slip  (open viscous wall flux)
  bool chi_wall_flux;              // chi_ac == dirichlet (open χ wall flux in divPol)
  bool core_flux_surface_average;  // firstX() && periodicY(xstart): automatic
};
```

Build it once in the constructor, make it a member, and print it with
`output_info.write`. With periodic X, `walls = false`, and none of the
radial logic runs.

### 3.3 Mode table

Applies when `walls` is true.

| Setting | `fixed` | `free_slip` | `no_slip` |
|---|---|---|---|
| Mean wall potential φ_b (stored in `Vort` radial guards) | dφ_b/dt = 0 | dφ_b/dt = 0 | Plan A law, κ = `phi_boundary_velocity` |
| φ̃ at wall | 0 | 0 | 0 |
| ψ̃ at wall; mean ψ, Pi_sum | 0; Neumann | 0; Neumann | 0; Neumann |
| Physical vorticity X guards (`vorticity`, used by stencils and viscosity) | Neumann | ω = 0 at midpoint | Neumann |
| Viscous vorticity and energy wall flux | closed | open | closed |
| Reynolds wall fluxes (`DivJReynolds`, `DivEReynolds`, Reynolds parts of `divPol`) | closed | closed | closed |
| E×B energy wall flux (`DivEExB`) | closed | closed | closed |
| `phiSolver` X flags | `INVERT_SET` | `INVERT_SET` | `INVERT_SET` |
| χ DC midpoint value | 0 | 0 | dφ_b/dt (Plan A §3.4) |
| χ AC | from `polarisation_ac_boundary` | same | same |
| χ wall flux in `divPol` | open only for `dirichlet` | same | same |
| Core boundary (closed surface) | flux-surface-averaged stored value | same | same, and flux-surface-averaged rate |

Physical meaning of each mode:

- **`fixed`:** the E×B potential at the wall is held. The slice's poloidal
  momentum is held by the boundary condition apart from changes in ψ̄_b. There
  is no viscous charge flux through the wall.
- **`free_slip`:** as `fixed`, plus the wall is stress-free: ω = 0 at the wall,
  and vorticity may diffuse viscously through it. In the stream
  function–vorticity form this is the free-slip condition.
- **`no_slip`:** the mean E×B flow at the wall relaxes to zero at velocity κ.
  This is a deliberate, resolution-independent torque.

`polarisation_ac_boundary` controls χ:

- **`neumann` (default):** ∂ₓχ̃ = 0 and the polarisation transport is closed.
  The boundary cells then have a small charge mismatch between ion density
  and vorticity, equal to the dropped wall flux. This mode allows larger
  timesteps.
- **`dirichlet`:** χ̃ = 0, consistent with Φ̃_wall = 0. The χ wall flux is
  open, so the polarisation drift carries ions and energy into the wall
  consistently with the vorticity.

The `neumann` condition is not fully insulating when g13 ≠ 0, because the
normal flux also contains g^{xz}∂_zχ̃. Document this; it is zero on
shifted-metric orthogonal grids.

### 3.4 Changes by function

1. **Constructor**
   - Read the enums (and `phi_boundary_velocity` only for `no_slip`), and
     build `RadialBoundaryModel`.
   - If `walls`:
     - Force `Vort:evolve_bndry` (throw if the user set it false).
     - Perform the `FieldData` swap and `initialiseBoundaryState()`.
     - Set `phiSolver` inner and outer flags to `INVERT_SET`.
     - Configure `chiSolver` flags: `INVERT_AC_GRAD` for `neumann`, 0 for
       `dirichlet` (both give DC midpoint Dirichlet). The DC value enters via
       the defect correction, not through the solver flags.
   - Throw if any of the user settings rejected in §2 are present
     (`isSet()`).
   - Keep the requirement `mesh->xstart >= 2` (two radial guard cells) for
     walls.
   - Remove the branches on `phi_boundary_relax`.
2. **`initialiseBoundaryState()`, `applyPotentialBoundary()`, `restartVars()`**
   - Make these unconditional when `walls` (no `phi_boundary_relax` checks).
   - Keep restart version 1 semantics.
   - Initial wall values are the midpoint of the initial φ (zero unless
     restarting). Document this.
3. **Wall rates (Plan A `computeWallRates`/`writeWallRates`)**
   - For `fixed` and `free_slip`, rates are zero and the χ DC value is zero.
   - For `no_slip`, use the Plan A law.
   - `finally()` always writes the radial storage cells of `ddt(Vort)` (zeros
     or rates) and writes interior cells from the physics.
4. **Pressure and ψ (`applyPressureBoundary`)**
   - Replace with one function, `applyWallPressureBoundary(Field3D& f)`,
     implementing the current average-z branch: ghost = 2·(z-mean of the
     first interior cell) − interior.
   - Fill **all** radial guard columns, not just two, by copying the first
     guard value outward.
   - Apply it to `Pi_sum` as now, **and** to `psi` after every
     `psi = pressure_weight * Pi_sum`. The weight's guard values come from
     species density boundaries, so ψ̃ = 0 is otherwise not guaranteed.
   - Delete the plain-Neumann branch and the option.
5. **Physical vorticity guards**
   - Add `applyVorticityWallBoundary(Field3D& w)`: X guards Neumann
     (`ghost = interior`) or, for `free_slip`, midpoint zero
     (`ghost = −interior`), filling all guard columns.
   - In `transform_impl`: `vorticity.applyBoundary()` (Y boundaries from
     user `Vort` settings), then `applyVorticityWallBoundary(vorticity)`,
     then communicate.
   - In `calculateViscosity()`, use `vorticity` directly; remove the
     hard-coded `"neumann"`.
6. **Wall fluxes**
   - Add a `bool` argument where needed and pass explicit values everywhere.
     Do not rely on defaults in this component.
   - Reynolds Q term in `DivJReynolds`:
     `Div_a_Grad_perp_XZ(0.5c, vEdotGradPsi, false)`.
   - `DivEReynolds`: `Div_a_Grad_perp_XZ(..., Phi, false)`, and rename the
     local `chi` to `Phi` with no boundary call (the physical guards of φ + ψ
     are already set).
   - Viscosity: `perpendicularFaceFluxes(nu, vorticity, i, model.viscous_wall_flux)`
     for both `DivJviscous` and `DivEviscous` (the same flag, so the charge
     and energy wall fluxes pair).
   - `divPol` χ term:
     `Div_a_Grad_perp_XZ(weight * c, chi, model.chi_wall_flux)`. All other
     `divPol` terms stay closed.
7. **Reynolds intermediates**
   - `vEdotGradPsi` (Q): X guards midpoint zero (ghost = −interior), which is
     exact because Q = 0 on an impermeable wall.
   - `DelpPhi_2B2`: X guards linear extrapolation (`free_o2`). These are used
     only by the reconstruction next to a closed face.
   - Y boundaries keep their current treatment.
8. **`applyPolarisationBoundary`**
   - Add an AC mode parameter: Neumann gives ghost AC = interior AC;
     Dirichlet gives ghost AC = −interior AC.
   - The DC midpoint value comes from the Plan A arrays.
   - Update the function's header comment.
9. **Core boundary**
   - Replace `phi_core_averagey` with `model.core_flux_surface_average`,
     computed automatically.
   - At closed flux surfaces the stored value and the rate are flux-surface
     averages (Plan A §3.3) in every mode.
10. **Header documentation**
    - Rewrite the option block to list only real options.
    - Fix the `poloidal_flows` default, and remove the stale entries listed
      in §1.2.

## 4. Implementation order

1. Enums, `RadialBoundaryModel`, option reading, rejection checks, and the
   startup summary. Compile.
2. Unify the wall-state path (remove the `phi_boundary_relax` branches).
   Rates are zero for `fixed` and `free_slip`. Compile and run regression
   §5.3.
3. Pressure and ψ wall boundary (§3.4 item 4).
4. Vorticity wall boundary and viscosity flags (§3.4 items 5 and 6).
5. Reynolds intermediates and fluxes (§3.4 items 6 and 7), and the
   shadowing rename.
6. χ AC option and `divPol` flux flag (§3.4 items 6 and 8).
7. Core flux-surface averaging cleanup (§3.4 item 9).
8. Documentation.

After each step: compile, run the unit tests, and run the §5.2 tests that
apply.

## 5. Tests

### 5.1 Unit tests (`tests/unit/`, GoogleTest, FakeMesh fixtures)

1. **Option handling.** Each enum value parses. An invalid string throws.
   `phi_boundary_velocity` is required for `no_slip`. For the other modes it
   is not read, so `options["phi_boundary_velocity"].valueUsed()` is false
   after construction, and BOUT++'s unused-option check rejects it at startup.
   Each rejected user setting in §2 throws.
2. **`applyWallPressureBoundary`.** The guard midpoint z-mean equals the
   interior z-mean (mean Neumann), and the guard's AC part is minus the
   interior AC (Dirichlet zero). All guard columns are filled.
3. **`applyVorticityWallBoundary`.** Neumann and midpoint-zero variants.
4. **`applyPolarisationBoundary`.** Both AC modes, with and without DC values.
5. **`perpendicularFaceFluxes` with `bndry_flux = false`.** The X-face flux at
   both walls is exactly zero, including the g13 corner term. Use a FakeMesh
   with non-zero g13.

### 5.2 Integrated tests

Use 2D x–z slabs, cloned from an existing integrated test's `runtest` and
`BOUT.inp`.

1. **Wall budgets (all three modes × both χ AC modes).** Run a seeded
   interchange-like or random-noise case with walls, closed parallel
   dynamics, no sources and zero viscosity (except a `free_slip` variant with
   viscosity). Each output reports:
   - **Ion particle content:** changes only through the χ wall flux for
     `dirichlet` (compute it from a diagnostic of that flux), and otherwise
     not at all.
   - **Total energy** (E_ion integral plus electron energy): closes to
     solver tolerance.
   - **Vorticity integral:** changes only through the viscous wall flux in
     `free_slip`.
2. **`no_slip` rigid rotation.** Repeat Plan A test 6.2 through
   `radial_boundary = no_slip`. τ_eff = L_phys/(2κ) at nx = 16, 32 and 64.
3. **`fixed` and `free_slip` hold the wall potential.** The stored φ_b is
   unchanged to round-off over a run.
4. **`free_slip` versus `fixed` with viscosity.** In a decaying shear-flow
   case, `free_slip` gives zero wall vorticity (diagnostic), and `fixed`
   gives zero normal gradient.
5. **Closed flux surfaces.** In a tokamak case with a core boundary, the
   applied core φ_b is identical across y to round-off in every mode.
6. **Restart round trip** in each mode.

### 5.3 Regression

Compare the old `phi_boundary_relax = true` with `neumann_boundary_average_z = true`
against the new `no_slip`, with κ = h_b/(2τ) chosen to match at one
resolution. Results should agree closely but not bitwise: the wall law is now
second order and the Reynolds and viscous wall fluxes are closed. Record the
differences in the PR description.

## 6. Known issues outside this plan

These should be raised with the maintainer and not changed silently.

- **No parallel current in the vorticity equation.**
  `ddt(Vort) = DivJReynolds + DivJdia + DivJviscous` has no
  `Div_par(Jpar_e + Jpar_i)` term and does not read `fields:DivJextra`. The
  energy equation does include the parallel electrical work (`DivEelectric`,
  `phi_Div_Je`). Hermes-3's `vorticity` component adds both. Unless this is
  deliberate for 2D tests, parallel currents change E_ion without changing
  ω, and sheath currents are absent from charge balance.
- **The Reynolds operator is still the original three-term split.** A
  skew-symmetric face-flux version with a corner stream function was
  prototyped separately. If adopted, its wall flux (keep the stream-function
  part at walls) and energy flux (K⊥ v_E) replace the corresponding
  entries in §3.3.
