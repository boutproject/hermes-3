# Polarisation drift and 3D energy: field contract

This contract describes field units, dimensions, numerical-face conventions
and ownership for the ion-energy and conservative polarisation components.
The equation and developer guides document the model and its regression
coverage. Geometry/analysis support limits below must be checked before
interpreting three-dimensional energy budgets.

## Units, dimensions and ownership

Let `P0 = qe Tnorm Nnorm`, `L0 = rho_s0`, `w0 = Omega_ci`, and
`v0 = Cs0`. Conversion multiplies stored normalised values to give SI.
Time converts with `1/w0`. Metrics and spacings must use their own conversion
attributes; their product `J dx dy dz` converts with `L0^3` for the slab
normalisation. Do not assign a Cartesian length conversion to tokamak angular
coordinates. Scalar normalisation fields and resolved settings are required.

`C` below means a cell-centred field with labelled dimensions `(t,x,y,z)`;
static metrics may omit time and broadcast over Z. xhermes may rename Y to
`theta`; adapters must discover its toroidal name. `Fx+` and `Fz+` mean arrays
indexed by the lower cell of the face: index `(i,j,k)` stores respectively
`(i+1/2,j,k)` and `(i,j,k+1/2)`. `Fx-` and `Fy-` are the BOUT++ numerical-flow
outputs at the lower face of the indexed cell. Opposite conventions require an
explicit index shift; names containing `xlow` must not be treated as `Fx+`.

All face currents/powers below already include numerical face area. Sum face
values directly, without multiplying by area or by `2*pi`. For a boundary
face oriented in the increasing coordinate direction, outward power has
negative sign at the lower boundary and positive sign at the upper boundary.
An inter-rank face, periodic connection or branch-cut connection is internal.
Only physical boundary owners contribute losses. Cell-centred source powers
are volume-integrated over physical cells, excluding all guards and topology
duplicates. Separate physical-target masks are required for multi-target
meshes; trimming the first and last Y rows is insufficient.

Integrals must declare the simulated sector and its angular extent. Full-torus
reporting multiplies every inventory, source and boundary exchange by the same
`full_torus_angle/sector_angle` factor, exactly once. Saved `dv` or integrated
face powers must first be checked against that domain convention.

## Existing reservoir and recovery output

| Field | Current owner | Meaning | SI units; conversion | Dimensions | Required use / boundary treatment |
| --- | --- | --- | --- | --- | --- |
| `Eion` | `evolve_vorticity_ion_energy` | Evolved ion thermal + parallel + perpendicular energy density | J m^-3; `P0` | C | Mandatory inventory; physical cells only |
| `Tion`, `P<s>` | Ion energy recovery (ion pressures in state), `evolve_pressure` for electrons | Common ion temperature and species pressure | eV; `Tnorm`, Pa; `P0` | C | Thermal reservoir `3/2 P`; ion pressure can be reconstructed from `N<s> Tion` |
| `N<s>` | `evolve_density` or dependent density component | Number density | m^-3; `Nnorm` | C | Mandatory with actual density owner recorded |
| `NV<s>` | `evolve_momentum` or flow closure | Mass-weighted parallel momentum | kg m^-2 s^-1; `Mp Nnorm v0` | C | Mandatory for each evolved parallel kinetic reservoir and boundary current |
| `V<s>` | Momentum recovery | Recovered parallel velocity | m s^-1; `v0` | C | Optional saved diagnostic; reconstruction must use the configured soft density floor |
| `ion_KE_parallel`, `ion_KE_perp` | Ion energy | Recovery-consistent kinetic densities | J m^-3; `P0` | C | Currently diagnose-only; saved perpendicular energy required for unsupported shifted-grid reconstruction |
| `phi` | Ion energy recovery, sheath modifies Y guards | Electrostatic potential | V; `Tnorm` | C | Completed-state potential needed for electrical work |
| `Vort` | Ion energy | Physical vorticity in output | C m^-3; `qe Nnorm` | C | Solver radial guards instead encode wall potential; do not integrate solver guard storage as vorticity |
| `ion_pressure_potential` | Ion energy, diagnose-only | `psi = sum A_s P_s/(Abar Z_s)` | V; `Tnorm` | C | Recovery and conservative transport input; pressure guards follow the wall model |
| `ion_chi`, `ddt(Phi)` | Conservative drift | Auxiliary current potential, approximate time derivative | V s^-1; `Tnorm w0` | C | `ion_chi` saved independently of diagnose; `ddt(Phi)` alias diagnose-only; selected radial AC / prescribed DC guards |
| `ion_phi_dissipation` | Ion finalisation, diagnose-only | Numerical parallel potential vorticity rate, excluded from chi and physical currents | A m^-3; `qe Nnorm w0` | C | Zero when disabled; physical cells only, preserving encoded radial rates; no explicit energy source |
| `ion_phi_iterations`, `ion_chi_iterations`, `ion_chi_residual` | Recovery for phi; conservative drift for chi | Recovery iteration count, fixed auxiliary solve count, max normalised FV residual | Dimensionless/count; residual conversion `qe Nnorm w0` if reported in A m^-3 | `(t)` | Implementation diagnostics; auxiliary residual must remain separate from missing power and time quadrature error |

Parallel kinetic density is `1/2 Mp A_s N_s V_s^2`, with the solver's
velocity/floor convention. Boussinesq perpendicular kinetic energy is counted
once for the common ion flow. Total plasma energy is `Eion + Ue + Kparallel_e`;
the independent sum of all ion/electron reservoirs is a consistency check.
Adding ion reservoirs to `Eion` would double count them.

## Stage-2 completed-current diagnostics

| Field | Owner | Meaning | SI units / conversion | Dimensions and use |
| --- | --- | --- | --- | --- |
| `parallel_current` | Conservative drift | Completed `sum (Z/A) momentum`, including sheath guards | A m^-2; `qe Nnorm v0` | C; current used by auxiliary RHS |
| `ion_parallel_current`, `electron_parallel_current` | Ion finalisation | Completed currents used by parallel electrical work | A m^-2; `qe Nnorm v0` | C; use with completed phi, including physical target guards |
| `DivJpar`, `ion_current_rhs`, `ion_auxiliary_residual` | Conservative drift | Parallel divergence, selected auxiliary RHS, RHS minus forward chi operator | A m^-3; `qe Nnorm w0` | C; physical cells including target-adjacent cells |
| `ddt(Vort)` | Ion finalisation | Completed vorticity rate including viscosity and optional numerical damping | A m^-3; `qe Nnorm w0` | C; physical cells only, exclude encoded radial solver guards |

These outputs require `diagnose=true`. The three divergence/RHS/residual fields
are also published in `fields:` independently of diagnose. The selected RHS is
`DivJReynolds+DivJdia+DivJpar+DivJextra+DivJcol`, omitting absent optional terms
and excluding viscosity and numerical damping. `DivJpol+DivJdia+
DivJpar+DivJextra+DivJcol` equals the selected RHS minus its forward chi operator.

Conservative drift reads completed charged momentum and phi with `Regions::All`,
establishing recovery -> production sheath -> drift -> finalisation. Periodic
X-Z planes reject a nonzero volume-weighted mean above their MPI-reduced
`128*epsilon*sum(abs(volume*rhs))` bound, without projecting it away.

In this BOUT revision, non-FCI `Field3D::setBoundaryTo` skips new physical Y
boundary regions. Restored solver `NV<s>` dumps consequently need not retain
completed sheath momentum guards. Use dedicated completed-current diagnostics
for target-current/electrical checks; zero solver momentum guards do not imply
zero target current. The unit regression checks species momenta directly in
the completed state.

## Existing transport and source output

| Field | Current owner | Meaning | SI units; conversion | Dimensions / orientation | Coverage and boundary ownership |
| --- | --- | --- | --- | --- | --- |
| `DivJReynolds`, `reynolds_work` | Ion energy | Divergence and face work of the retained Reynolds method | A m^-3; `qe Nnorm w0`, W m^-3; `P0 w0` | C | Diagnose-only output; Reynolds faces and divergence are always published in state |
| `DivJdia` | Ion energy state | Total diamagnetic-current divergence | A m^-3; `qe Nnorm w0` | C | Required conservative-current input; do not infer zero from an absent optional diagnostic |
| `polarisation_flux_x`, `polarisation_flux_z` | Conservative drift | `-JReynolds - c grad(chi)`, area-integrated X-Z current | A; `qe Nnorm w0 L0^3` | Fx+, Fz+ | Diagnose-only; physical X faces retained; optional poloidal extension excluded from these arrays |
| `DivJpol` | Conservative drift | Full conservative current divergence | A m^-3; `qe Nnorm w0` | C | Includes physical wall current and poloidal extension; legacy drift uses the same name for `-DivJ`, requiring implementation/version metadata |
| `DivEExB`, `DivEelectric`, `DivEparallel`, `DivEpol` | Ion energy; conservative drift owns `DivEpol` | Numerical energy-flux divergences used by Eion RHS | W m^-3; `P0 w0` | C | Diagnose-only; spatial integral gives total signed boundary exchange, but cannot separate individual targets/walls without faces |
| `ddt(Eion)` | Ion energy | Instantaneous evolved energy RHS | W m^-3; `P0 w0` | C | Diagnose-only; independent spatial/source-closure audit |
| `phi_Div_Je`, `electron_exb_work` | Ion energy | Electron electrical and diamagnetic work coupled to ions | W m^-3; `P0 w0` | C | Internal exchange only with the corresponding evolved electron reservoir |
| `zero_current_defect` (analysis API), `model_defect_power`, `known_model_defect_power` (analysis results) | Analysis of electron pressure / ion electrical work | Unpaired electron force work `Ve Grad_par(Pe) + Je Grad_par(phi)` | W; volume integral of normalised work times `P0 w0 L0^3` | `(t)` | Explicit closure defect; separate from heating. Uniform-Y metrics and C2 gradients with valid guards or declared periodicity; nonzero/unsupported work precludes a claimed closed balance |
| `DivJviscous`, `DivEviscous`, `kinematic_viscous_heating`, `electron_viscous_heating` | Ion energy | Kinematic current, diagnostic stress-energy flux and heating partition | A m^-3 or W m^-3; respective current/power conversion above | C | Viscous drift work already carried by `DivEpol`; `DivEviscous` is not an additional loss |
| `phi_wall_rate` | Ion energy | Encoded midpoint potential rate | V s^-1; `Tnorm w0` | C radial guards | Diagnose-only; only physical wall rows, zero elsewhere; core surface averaging follows connected-Y ownership |
| `polarisation_wall_*weight`, `polarisation_wall_pressure_mismatch` | Conservative drift | Uniform-Z transport work/thermal weights and pressure mismatch | V; `Tnorm` | Fx+ on physical X walls | Diagnose-only; use with the same face current, do not substitute an unaveraged cell product |
| `particle_pol_wall_loss_<s>` | Conservative drift | Outward polarisation particle loss per adjacent-cell volume | m^-3 s^-1; `Nnorm w0` | C | Diagnose-only; volume-integrate exactly once |
| `Eion_src` | Ion energy | Converted total energy sources from components | W m^-3; `P0 w0` | C | Contains multiple classifications; not directly equal to external input power |
| `S<s>_src` | Density evolution | Applied external density source, including mask/time prefactor | m^-3 s^-1; `Nnorm w0` | C | Diagnose-only; preferable to reconstructing a copied analytic source shape |
| `P<s>_src` | Pressure evolution | Applied pressure source | Pa s^-1; `P0 w0` | C | Diagnose-only; power density is `3/2` this field; removed ion pressure evolution removes its source owner |
| `SN<s>`, `SNV<s>`, `SP<s>` | Density, momentum, pressure evolution | Aggregated component sources | Corresponding density/momentum/pressure rate | C | Mixed ownership; retain individual external/internal/sheath contribution fields before aggregating |
| `pf<s>_tot_xlow/ylow`, `mf<s>_tot_xlow/ylow` | Density / momentum evolution | Numerical particle/momentum face transport | s^-1; `Nnorm w0 L0^3`, N; `Mp Nnorm v0 w0 L0^3` | Fx-, Fy- | Diagnose-only and allocation-dependent; momentum transport alone is not kinetic-energy power |
| `ef<s>_tot_xlow/ylow`, `ef<s>_adv_ylow`, `ef<s>_visc_heat_ylow` | Pressure evolution | Numerical thermal-energy faces, advection and numerical viscous heating flux | W; `P0 w0 L0^3` | Fx-, Fy- | Explicitly may be incomplete; inspect selected components and sum each contribution once |
| `E<species-pair>_coll` | `braginskii_heat_exchange` | Pairwise thermal exchange | W m^-3; `P0 w0` | C | Diagnose-only; volume-integrated pair contributions must cancel |

Version 5 of `evolve_vorticity_ion_energy` pairs with version 1 of
`polarisation_drift`, with `implementation = conservative`. The auxiliary
options `polarisation_ac_boundary`, `polarisation_corrections` and `laplacian`
are owned by the drift section. Version-4 dumps retain their original
ion-energy option ownership. Version 6 pairs with conservative drift version 2,
adds completed parallel/optional currents to the auxiliary and vorticity RHS,
and evaluates parallel electrical work after sheath conditioning. Version 7
adds optional collective-sound-speed potential damping during ion finalisation;
the auxiliary equation is unchanged. Analysis supports versions 4/5/6/7. Standalone
`polarisation_drift` defaults to `legacy`.

The following state fields are published on every RHS independently of
`diagnose`; they are component inputs, with normalised values and valid
physical cells plus the guards/faces needed by the existing X-Z operators:

| State field | Owner | Meaning and dimensions |
| --- | --- | --- |
| `fields:ion_boussinesq_coefficient` | Ion recovery | Static metric coefficient `Abar/B²`; includes guards |
| `fields:ion_pressure_potential`, `fields:ion_mass_density` | Ion recovery | Cell pressure potential and actual ion mass density |
| `fields:ion_reynolds_x`, `fields:ion_reynolds_z` | Ion recovery | Area-integrated increasing-coordinate Reynolds faces, including owned X-wall faces |
| `fields:DivJReynolds`, `fields:DivJdia` | Ion recovery | Cell current divergences |
| `fields:ion_wall_rate` | Ion recovery | Encoded DC wall-potential rates in radial guards; zero elsewhere |
| `fields:ion_polarisation_energy_divergence`, `fields:DivJpol`, `fields:ion_chi` | Conservative drift | Cell energy/current divergences and auxiliary potential |
| `species:<s>:polarisation_density_source`, `species:<s>:polarisation_momentum_source` | Conservative drift | Separately identified particle/momentum transport sources; momentum field exists when species momentum exists |

Ion finalisation subtracts the kinetic conversion of these identified
polarisation sources because their parallel kinetic transport is already
included in `DivEpol`. This replaces the former private source cache without
changing the arithmetic. Source fields are not external heating.

`sheath_boundary` adds heat corrections to species energy sources in adjacent
physical Y cells, after conditioning momentum, velocity and potential guards.
It currently does not retain separate correction outputs. Its correction must
be classified as outward target power, using minus the volume integral of the
added source. A gamma-based total sheath power cannot be added on top of
numerical advective power and the additional correction.

## Required additions before a claimed complete 3D budget

| Required quantity | Proposed owner / storage | Integration and classification |
| --- | --- | --- |
| Boussinesq coefficient, psi, Reynolds X-Z faces/divergence, diamagnetic divergence and DC wall rates | Ion energy -> state, independent of `diagnose` | Component inputs; keep Reynolds method and evaluation in its existing owner |
| Auxiliary solution, current residual including all enabled physical currents, `DivEpol`, polarisation species source contributions | Conservative drift -> state; retained outputs | Current identity checked in physical boundary cells and interior; energy source cache replaced by state fields |
| `sum (Z_s/A_s) momentum_s`, `DivJpar`, enabled `DivJextra`, `DivJcol` | Completed-state current reads after sheath boundary | Include identical physical terms in auxiliary RHS and vorticity RHS; auxiliary excludes viscosity and numerical damping |
| Applied ion heating independent of ion pressure evolution | Explicit heating component -> applied SI source and integrated power | External input; core mask and time prefactor recorded; MMS forcing excluded |
| Numerical ion energy X/Y faces, separated by ExB, polarisation, enthalpy/kinetic and electrical transport | Respective flux owner -> retained area-integrated face powers or owner-tagged boundary histories | Signed power for each physical X wall and each Y target; electrical terms use completed sheath potential |
| Electron kinetic transport and pressure work, electrical exchange and frictional mechanical work | Electron momentum/pressure and exchange owners | Require evolved receiving reservoir; check internal sum, retain explicit model/numerical defects |
| Sheath correction for each species and target | Sheath -> separate added-source diagnostic or integrated target history | Boundary exchange; count once alongside numerical advective power |
| Enabled conduction/diffusion/radiation/sinks and wall work | Actual numerical term owner | Saved numerical face or source histories; closure follows what the RHS uses, including conduction suppression where sheath heat is imposed |
| Physical cell/target masks and sector volume | xhermes topology/geometry with minimal labelled 3D adapter | Cross-check solver boundary ownership; reject unsupported geometry rather than using a Cartesian perpendicular gradient |
| Saved recovery-consistent perpendicular energy | Ion energy | Mandatory for shifted grids with unsupported independent face reconstruction; provenance explicit |

Every enabled unmeasured term leaves coverage incomplete. Keep known histories,
list the missing term and omit a closed residual. With complete coverage,
positive outward power gives `rP = dE/dt - Pin + PXout + PYout + Pother_loss`.
The primary transient test is the cumulative energy change minus the integral
of signed net input power. Numerical defects are separate ledger entries.
Sampling error, integrated instantaneous-RHS error and auxiliary-current
residual are separate measurements.

## xhermes API audit

The configured `.venv/bin/python` resolves xhermes through its installed
editable package; no source-path injection is needed. `xhermes.open` and
`open_hermesdataset` exist in `xhermes/load.py`. The accessor has radial/poloidal
selectors and two-/four-target topology metadata. Its tokamak geometry method
assigns `dv`, `dr`, `hthe` and other fields with explicit `['x', 'theta']`
dimensions and eager `.data` products. That path needs a tested labelled 3D
adapter before it can handle a toroidal metric dimension safely. The local
energy-analysis loader currently rejects non-identity transforms and tokamak
topology. None of these gaps justifies copying a second topology reader or
editing `external/xhermes`.

## Validation limits

Integrated/MMS criteria declare the effective acceptance limits. Keep
recovery/auxiliary residuals, physical identities, spatial error, temporal
comparison and MPI comparison separate. Chi magnitude is diagnostic-only in
the parallel-current MMS, while its spatial orders and temporal/MPI checks
remain required. Numerical damping is excluded from auxiliary continuity;
full vorticity-rate checks subtract damping, viscosity and any injected MMS
forcing. Component-reordering trajectory comparisons are omitted after their
near-roundoff differences were accepted. The equation/developer guides explain
the distinct zero-curvature and curved-wall Reynolds convergence expectations.
