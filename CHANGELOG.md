# Changelog

All notable changes to Hermes-3 are documented in this file. Versions correspond
to git tags; the docker images (`boutproject/hermes-3`, `-builder`, `-jupyter`)
are built when the matching GitHub Release is published.

## v1.5.0

Changes since v1.4.1 (tagged 2025-12-09). ~122 pull requests merged.

### Physics & model features

- 1D fieldline geometry: reformulated cross-field transport with new geometry
  terms and diagnostics (#232)
- Sheath boundary features, including supersonic handling; `sin_alpha` as a
  `Field3D` (#574, #377)
- External Apar flutter for electromagnetic runs (#554)
- Target pump (#491)
- Viscous heating and energy conservation (#455)
- `neutral_full_velocity` transport coefficient (#454); optional nonorthogonal
  operators in `neutral_mixed` (#497)
- Numerical method choices (#578); relax-potential component (#465)

### Reactions & component framework

- Metadata attached to component objects (#596)
- State-variable access control in `transform()` methods (#421)
- Decouple reaction data (#534); reaction classes added to the `hermes`
  namespace (#538)
- Closure refactored to allow other forms (#399); H-isotope charge exchange
  moved to the new reactions framework (#398)

### Bug fixes

- Guard-cell communication fixes in `vorticity` and `braginskii_ion_viscosity`
  (#589)
- Fast-recycling heat flux (#494); `relax_potential` radial boundary condition
  (#528); `sheath_boundary_simple` (#482)
- Numerous option-permission / `AA` fixes (#540, #517, #509, #468, #466, #446)

### BOUT++ & build

- Adopt BOUT++ template (lazy) expressions (#585)
- BOUT++ version bumps (#440, #489, #563)
- Spack builds BOUT++ as a separate package (#376, #479, #488); spack builds of
  hermes-3 / BOUT-dev default to `RelWithDebInfo`, and BOUT-Spack updated to fix
  python dependency versions for tests (#619)
- Minimum CMake version set to 3.25 (#514)

### Testing

- More unit tests (#483); energy conservation tests (#455)
- 2D integrated CI test (#473); D-T-He 1D-recycling regression test (#425);
  nonorthogonal MMS test (#285)
- Faster and cleaned-up tests (#533, #610, #611, #612)

### Documentation & CI/infrastructure

- Docs: `diamagnetic_drift` description (#558), installation instructions
  (#569, #505), sphinx / `neutral_mixed` equation fixes (#516, #458)
- CI/docker: multi-arch images (#601), auto-delete old images (#606), stale-bot
  (#607, #617), tidy actions (#616), clang-format / python-format checks
  (#485, #486, #510, #545), automatic Python-dependency updates (#582, #618), PR
  template (#595)
