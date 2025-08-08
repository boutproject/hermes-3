#include <algorithm>

#include "../include/hermes_utils.hxx"
#include "../include/sheath_boundary.hxx"

#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/field.hxx>
#include <bout/field3d.hxx>
#include <bout/mesh.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/region.hxx>

#include <fmt/core.h>

#if __cpp_lib_unreachable >= 202202L
// C++23
#include <utility>
#define HERMES_UNREACHABLE std::unreachable()
#else
#include <stdexcept>
#define HERMES_UNREACHABLE throw std::logic_error("unreachable")
#endif

using bout::globals::mesh;

namespace {
/// Limited free gradient of log of a quantity
/// This ensures that the guard cell values remain positive
/// while also ensuring that the quantity never increases
///
///  fm  fc | fp
///         ^ boundary
///
/// exp( 2*log(fc) - log(fm) )
BoutReal limitFree(BoutReal fm, BoutReal fc, SheathLimitMode mode) {
  if ((fm < fc) && (mode == SheathLimitMode::limit_free)) {
    return fc; // Neumann rather than increasing into boundary
  }
  if (fm < 1e-10) {
    return fc; // Low / no density condition
  }

  BoutReal fp = 0;
  switch (mode) {
  case SheathLimitMode::limit_free:
  case SheathLimitMode::exponential_free:
    fp = SQ(fc) / fm; // Exponential
  case SheathLimitMode::linear_free:
    fp = 2.0 * fc - fm; // Linear
  }

#if CHECKLEVEL >= 2
  if (!std::isfinite(fp)) {
    throw BoutException("SheathBoundary limitFree: {}, {} -> {}", fm, fc, fp);
  }
#endif

  return fp;
}

/// Strong type to refer to a given boundary
enum class YBoundarySide : bool {
  lower, // ystart boundary
  upper, // yend boundary
};

/// Information about a point on a Y boundary
struct YBoundary {
  YBoundarySide sense; //< Is this boundary lower or upper?
  Ind3D i;             //< 3D index at this point
  Ind3D ip;            //< 3D index one step into the boundary
  Ind3D im;            //< 3D index one step away from the boundary
  int toward{};        //< Which direction in y is toward the boundary?
  YBoundary(Ind3D i, YBoundarySide sense)
      : sense(sense), i(i), ip(is_lower() ? i.yp() : i.ym()),
        im(is_lower() ? i.ym() : i.yp()), toward(is_lower() ? -1 : +1) {}
  /// Returns true if this point is on the lower Y boundary
  bool is_lower() const { return sense == YBoundarySide::lower; }
};

/// Loop over the ``sense`` Y boundary, calling ``func`` at each point
/// with a `YBoundary`
template <YBoundarySide sense, class Func>
void loop_boundary(const Mesh& mesh, Func func) {
  const bool is_lower = (sense == YBoundarySide::lower);
  auto range = is_lower ? mesh.iterateBndryLowerY() : mesh.iterateBndryUpperY();

  for (; !range.isDone(); ++range) {
    for (int jz = 0; jz < mesh.LocalNz; jz++) {
      // Y-index at this boundary
      const auto y = is_lower ? mesh.ystart : mesh.yend;
      // 3D index at this point
      const auto ny = mesh.LocalNy;
      const auto nz = mesh.LocalNz;
      const auto i = Ind3D{((range.ind * ny + y) * nz) + jz, ny, nz};
      const YBoundary bdry{i, sense};
      func(bdry);
    }
  }
}

/// Calculate potential phi assuming zero current
Field3D calculate_phi_simple(Options& allspecies, bool lower_y, bool upper_y,
                             const Field3D& Ne, const Field3D& Te, BoutReal Me,
                             BoutReal Ge, const Field3D& wall_potential,
                             SheathLimitMode density_boundary_mode,
                             SheathLimitMode temperature_boundary_mode,
                             BoutReal sheath_ion_polytropic) {

  // Need to sum  n_i Z_i C_i over all ion species
  //
  // To avoid looking up species for every grid point, this
  // loops over the boundaries once per species.
  Field3D ion_sum = 0.0;

  // Iterate through charged ion species
  for (const auto& kv : allspecies.getChildren()) {
    const Options& species = kv.second;

    if ((kv.first == "e") or !species.isSet("charge")
        or (get<BoutReal>(species["charge"]) == 0.0)) {
      continue; // Skip electrons and non-charged ions
    }

    const Field3D Ni = getNoBoundary<Field3D>(species["density"]);
    const Field3D Ti = getNoBoundary<Field3D>(species["temperature"]);
    const BoutReal Mi = getNoBoundary<BoutReal>(species["AA"]);
    const BoutReal Zi = getNoBoundary<BoutReal>(species["charge"]);
    Field3D Vi = species.isSet("velocity")
                     ? toFieldAligned(getNoBoundary<Field3D>(species["velocity"]))
                     : zeroFrom(Ni);

    const auto set_ion_sum_boundary = [&](const YBoundary& bdry) {
      const auto i = bdry.i;
      const auto ip = bdry.ip;
      // Free gradient of log density and temperature
      // This ensures that the guard cell values remain positive
      const BoutReal Ni_im = limitFree(Ni[ip], Ni[i], density_boundary_mode);
      const BoutReal Ti_im = limitFree(Ti[ip], Ti[i], temperature_boundary_mode);
      const BoutReal Te_im = limitFree(Te[ip], Te[i], temperature_boundary_mode);

      // Calculate sheath values at half-way points (cell edge)
      const BoutReal nisheath = 0.5 * (Ni_im + Ni[i]);
      // electron temperature
      const BoutReal tesheath = floor(0.5 * (Te_im + Te[i]), 1e-5);
      // ion temperature
      const BoutReal tisheath = floor(0.5 * (Ti_im + Ti[i]), 1e-5);

      // Sound speed
      const BoutReal C_i = sqrt((sheath_ion_polytropic * tisheath + Zi * tesheath) / Mi);
      const BoutReal visheath =
          bdry.is_lower() ? std::min(Vi[i], -C_i) : std::max(Vi[i], C_i);

      ion_sum[i] += bdry.toward * Zi * nisheath * visheath;
    };

    if (lower_y) {
      loop_boundary<YBoundarySide::lower>(*Ni.getMesh(), set_ion_sum_boundary);
    }

    if (upper_y) {
      loop_boundary<YBoundarySide::upper>(*Ni.getMesh(), set_ion_sum_boundary);
    }
  }

  Field3D phi = emptyFrom(Ne);

  const auto set_phi_boundary = [&](const YBoundary& bdry) {
    const auto i = bdry.i;
    const auto ip = bdry.ip;

    const BoutReal Ne_im = limitFree(Ne[ip], Ne[i], density_boundary_mode);
    const BoutReal Te_im = limitFree(Te[ip], Te[i], temperature_boundary_mode);

    // Calculate sheath values at half-way points (cell edge)
    const BoutReal nesheath = 0.5 * (Ne_im + Ne[i]);
    const BoutReal tesheath = floor(0.5 * (Te_im + Te[i]), 1e-5);

    phi[i] =
        tesheath * log(sqrt(tesheath / (Me * TWOPI)) * (1. - Ge) * nesheath / ion_sum[i]);

    // Add bias potential
    phi[i] += wall_potential[i];
    // Constant into sheath
    phi[i.yp()] = phi[i.ym()] = phi[i];
  };

  // ion_sum now contains the ion current, sum Z_i n_i C_i over all ion species
  // at mesh->ystart and mesh->yend indices
  if (lower_y) {
    loop_boundary<YBoundarySide::lower>(*phi.getMesh(), set_phi_boundary);
  }

  if (upper_y) {
    loop_boundary<YBoundarySide::upper>(*phi.getMesh(), set_phi_boundary);
  }

  return phi;
}

/// Calculate potential phi assuming zero current.
/// Note: This is equation (22) in Tskhakaya 2005, with I = 0
Field3D calculate_phi_Tskhakaya(Options& allspecies, bool lower_y, bool upper_y,
                                const Field3D& Ne, const Field3D& Te, BoutReal Me,
                                BoutReal Ge, const Field3D& wall_potential,
                                BoutReal sin_alpha) {

  // Need to sum  s_i Z_i C_i over all ion species
  //
  // To avoid looking up species for every grid point, this
  // loops over the boundaries once per species.
  Field3D ion_sum{zeroFrom(Ne)};

  // Iterate through charged ion species
  for (const auto& kv : allspecies.getChildren()) {
    Options& species = allspecies[kv.first];

    if ((kv.first == "e") or !IS_SET(species["charge"])
        or (get<BoutReal>(species["charge"]) == 0.0)) {
      continue; // Skip electrons and non-charged ions
    }

    const Field3D Ni =
        toFieldAligned(floor(GET_NOBOUNDARY(Field3D, species["density"]), 0.0));
    const Field3D Ti = toFieldAligned(GET_NOBOUNDARY(Field3D, species["temperature"]));
    const BoutReal Mi = GET_NOBOUNDARY(BoutReal, species["AA"]);
    const BoutReal Zi = GET_NOBOUNDARY(BoutReal, species["charge"]);

    const BoutReal adiabatic = IS_SET(species["adiabatic"])
                                   ? get<BoutReal>(species["adiabatic"])
                                   : 5. / 3; // Ratio of specific heats (ideal gas)

    const auto set_ion_sum_boundary = [&](const YBoundary& bdry) {
      const auto i = bdry.i;
      const auto ip = bdry.ip;

      // Free boundary extrapolate ion concentration
      // Limit range to [0,1]
      BoutReal s_i = std::clamp(0.5 * (3. * Ni[i] / Ne[i] - Ni[ip] / Ne[ip]), 0.0, 1.0);

      if (!std::isfinite(s_i)) {
        s_i = 1.0;
      }
      const BoutReal te = Te[i];
      const BoutReal ti = Ti[i];

      // Equation (9) in Tskhakaya 2005
      BoutReal grad_ne = -bdry.toward * (Ne[ip] - Ne[i]);
      BoutReal grad_ni = -bdry.toward * (Ni[ip] - Ni[i]);

      // Note: Needed to get past initial conditions, perhaps transients
      // but this shouldn't happen in steady state
      // Note: Factor of 2 to be consistent with later calculation
      if (std::abs(grad_ni) < 2e-3) {
        grad_ni = grad_ne = 2e-3; // Remove kinetic correction term
      }

      // Limit for e.g. Ni zero gradient
      const BoutReal C_i_sq =
          std::clamp((adiabatic * ti + Zi * s_i * te * grad_ne / grad_ni) / Mi, 0., 100.);

      // Note: Vzi = C_i * sin(α)
      ion_sum[i] += s_i * Zi * sin_alpha * sqrt(C_i_sq);
    };

    if (lower_y) {
      loop_boundary<YBoundarySide::lower>(*mesh, set_ion_sum_boundary);
    }

    if (upper_y) {
      loop_boundary<YBoundarySide::upper>(*mesh, set_ion_sum_boundary);
    }
  }

  Field3D phi = emptyFrom(Ne); // So phi is field aligned

  const auto set_phi_boundary = [&](const YBoundary& bdry) {
    const auto i = bdry.i;
    if (Te[i] <= 0.0) {
      phi[i] = 0.0;
    } else {
      phi[i] = Te[i] * log(sqrt(Te[i] / (Me * TWOPI)) * (1. - Ge) / ion_sum[i]);
    }

    phi[i] += wall_potential[i];        // Add bias potential
    phi[i.yp()] = phi[i.ym()] = phi[i]; // Constant into sheath
  };

  // ion_sum now contains  sum  s_i Z_i C_i over all ion species
  // at mesh->ystart and mesh->yend indices
  if (lower_y) {
    loop_boundary<YBoundarySide::lower>(*phi.getMesh(), set_phi_boundary);
  }

  if (upper_y) {
    loop_boundary<YBoundarySide::upper>(*phi.getMesh(), set_phi_boundary);
  }
  return phi;
}

/// If phi is set, use free boundary condition;
/// If phi not set, calculate assuming zero current, according to sheath kind
template <hermes::SheathKind kind>
Field3D init_phi(Options& state, Options& allspecies, bool lower_y, bool upper_y,
                 const Field3D& Ne, const Field3D& Te, BoutReal Me, BoutReal Ge,
                 const Field3D& wall_potential, BoutReal sin_alpha,
                 SheathLimitMode density_boundary_mode,
                 SheathLimitMode temperature_boundary_mode,
                 BoutReal sheath_ion_polytropic) {

  if (IS_SET_NOBOUNDARY(state["fields"]["phi"])) {
    return toFieldAligned(getNoBoundary<Field3D>(state["fields"]["phi"]));
  }

  if constexpr (kind == hermes::SheathKind::normal) {
    return calculate_phi_Tskhakaya(allspecies, lower_y, upper_y, Ne, Te, Me, Ge,
                                   wall_potential, sin_alpha);
  } else if constexpr (kind == hermes::SheathKind::simple) {
    return calculate_phi_simple(allspecies, lower_y, upper_y, Ne, Te, Me, Ge,
                                wall_potential, density_boundary_mode,
                                temperature_boundary_mode, sheath_ion_polytropic);
  } else if constexpr (kind == hermes::SheathKind::insulating) {
    // No model for potential for this case, but return something to
    // avoid duplicated checks later
    return 0.0;
  } else {
    static_assert(false, "Unhandled sheath kind");
  }
}

} // namespace

template <hermes::SheathKind kind>
SheathBoundaryBase<kind>::SheathBoundaryBase(Options& options, BoutReal Tnorm) {
  AUTO_TRACE();

  Ge = options["secondary_electron_coef"]
           .doc("Effective secondary electron emission coefficient")
           .withDefault(0.0);

  if ((Ge < 0.0) or (Ge > 1.0)) {
    throw BoutException("Secondary electron emission must be between 0 and 1 ({:e})", Ge);
  }

  sin_alpha = options["sin_alpha"]
                  .doc("Sin of the angle between magnetic field line and wall surface. "
                       "Should be between 0 and 1")
                  .withDefault(1.0);

  if ((sin_alpha < 0.0) or (sin_alpha > 1.0)) {
    throw BoutException("Range of sin_alpha must be between 0 and 1");
  }

  gamma_e = options["gamma_e"]
                .doc("Electron sheath heat transmission coefficient")
                .withDefault(3.5);

  gamma_i =
      options["gamma_i"].doc("Ion sheath heat transmission coefficient").withDefault(3.5);

  sheath_ion_polytropic = options["sheath_ion_polytropic"]
                              .doc("Ion polytropic coefficient in Bohm sound speed")
                              .withDefault(1.0);

  lower_y = options["lower_y"].doc("Boundary on lower y?").withDefault<bool>(true);
  upper_y = options["upper_y"].doc("Boundary on upper y?").withDefault<bool>(true);

  always_set_phi =
      options["always_set_phi"]
          .doc("Always set phi field? Default is to only modify if already set")
          .withDefault<bool>(false);

  // Read wall voltage, convert to normalised units
  wall_potential = options["wall_potential"]
                       .doc("Voltage of the wall [Volts]")
                       .withDefault(Field3D(0.0))
                   / Tnorm;
  // Convert to field aligned coordinates
  wall_potential = toFieldAligned(wall_potential);
  // Note: wall potential at the last cell before the boundary is used,
  // not the value at the boundary half-way between cells. This is due
  // to how twist-shift boundary conditions and non-aligned inputs are
  // treated; using the cell boundary gives incorrect results.

  floor_potential = options["floor_potential"]
                        .doc("Apply a floor to wall potential when calculating Ve?")
                        .withDefault<bool>(true);

  no_flow = options["no_flow"]
                .doc("Set zero particle flow, keeping energy flow")
                .withDefault<bool>(false);

  density_boundary_mode = options["density_boundary_mode"]
                              .doc("BC mode: limit_free, exponential_free, linear_free")
                              .withDefault<SheathLimitMode>(SheathLimitMode::limit_free);

  pressure_boundary_mode = options["pressure_boundary_mode"]
                               .doc("BC mode: limit_free, exponential_free, linear_free")
                               .withDefault<SheathLimitMode>(SheathLimitMode::limit_free);

  temperature_boundary_mode =
      options["temperature_boundary_mode"]
          .doc("BC mode: limit_free, exponential_free, linear_free")
          .withDefault<SheathLimitMode>(SheathLimitMode::limit_free);

  diagnose = options["diagnose"]
                 .doc("Save additional output diagnostics")
                 .withDefault<bool>(false);
}

template <hermes::SheathKind kind>
void SheathBoundaryBase<kind>::transform(Options& state) {
  AUTO_TRACE();

  Options& allspecies = state["species"];
  Options& electrons = allspecies["e"];

  // Need electron properties
  // Not const because boundary conditions will be set
  Field3D Ne = toFieldAligned(floor(GET_NOBOUNDARY(Field3D, electrons["density"]), 0.0));
  Field3D Te = toFieldAligned(GET_NOBOUNDARY(Field3D, electrons["temperature"]));
  Field3D Pe = IS_SET_NOBOUNDARY(electrons["pressure"])
                   ? toFieldAligned(getNoBoundary<Field3D>(electrons["pressure"]))
                   : Te * Ne;

  // Ratio of specific heats
  const BoutReal electron_adiabatic =
      (kind == hermes::SheathKind::normal and IS_SET(electrons["adiabatic"]))
          ? get<BoutReal>(electrons["adiabatic"])
          : 5. / 3;

  // Mass, normalised to proton mass
  const BoutReal Me =
      IS_SET(electrons["AA"]) ? get<BoutReal>(electrons["AA"]) : SI::Me / SI::Mp;

  // This is for applying boundary conditions
  Field3D Ve = IS_SET_NOBOUNDARY(electrons["velocity"])
                   ? toFieldAligned(getNoBoundary<Field3D>(electrons["velocity"]))
                   : zeroFrom(Ne);

  Field3D NVe = IS_SET_NOBOUNDARY(electrons["momentum"])
                    ? toFieldAligned(getNoBoundary<Field3D>(electrons["momentum"]))
                    : zeroFrom(Ne);

  Coordinates* coord = mesh->getCoordinates();

  //////////////////////////////////////////////////////////////////
  // Electrostatic potential
  Field3D phi = init_phi<kind>(state, allspecies, lower_y, upper_y, Ne, Te, Me, Ge,
                               wall_potential, sin_alpha, density_boundary_mode,
                               temperature_boundary_mode, sheath_ion_polytropic);

  // Field to capture total sheath heat flux for diagnostics
  Field3D electron_sheath_power_ylow = zeroFrom(Ne);

  //////////////////////////////////////////////////////////////////
  // Electrons

  Field3D electron_energy_source =
      electrons.isSet("energy_source")
          ? toFieldAligned(getNonFinal<Field3D>(electrons["energy_source"]))
          : zeroFrom(Ne);

  hflux_e = zeroFrom(electron_energy_source); // sheath heat flux for diagnostics

  // Cross-sectional area in XZ plane and cell volume.
  // Would be nice to cache these
  auto da_m = emptyFrom(coord->J);    // [m^2]
  auto da_p = emptyFrom(coord->J);    // [m^2]
  auto dv = emptyFrom(coord->J);      // [m^3]
  auto da_dv_m = emptyFrom(coord->J); // [m^-1]
  auto da_dv_p = emptyFrom(coord->J); // [m^-1]
  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {
      auto i = indexAt(Ne, r.ind, mesh->ystart, jz);
      auto im = i.ym();
      da_m[i] = (coord->J[i] + coord->J[im])
                / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[im])) * 0.5
                * (coord->dx[i] + coord->dx[im]) * 0.5 * (coord->dz[i] + coord->dz[im]);

      dv[i] = (coord->dx[i] * coord->dy[i] * coord->dz[i] * coord->J[i]);

      da_dv_m[i] = da_m[i] / dv[i];
    }
  }

  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {
      auto i = indexAt(Ne, r.ind, mesh->yend, jz);
      auto ip = i.yp();
      da_p[i] = (coord->J[i] + coord->J[ip])
                / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[ip])) * 0.5
                * (coord->dx[i] + coord->dx[ip]) * 0.5
                * (coord->dz[i] + coord->dz[ip]); // [m^2]

      dv[i] = (coord->dx[i] * coord->dy[i] * coord->dz[i] * coord->J[i]);

      da_dv_p[i] = da_p[i] / dv[i];
    }
  }

  const auto set_electron_boundaries = [&](const YBoundary& bdry) {
    const auto i = bdry.i;
    const auto im = bdry.im;
    const auto ip = bdry.ip;

    // Free gradient of log electron density and temperature
    // Limited so that the values don't increase into the sheath
    // This ensures that the guard cell values remain positive
    Ne[im] = limitFree(Ne[ip], Ne[i], density_boundary_mode);
    Te[im] = limitFree(Te[ip], Te[i], temperature_boundary_mode);
    Pe[im] = limitFree(Pe[ip], Pe[i], pressure_boundary_mode);

    // Free boundary potential linearly extrapolated
    phi[im] = 2 * phi[i] - phi[ip];

    if constexpr (is_insulating()) {
      // Set zero flow through boundary
      // This will be modified when iterating over the ions
      Ve[bdry.im] = -Ve[bdry.i];
      NVe[bdry.im] = -NVe[bdry.i];
      return;
    }

    const BoutReal nesheath = 0.5 * (Ne[im] + Ne[i]);
    const BoutReal tesheath = 0.5 * (Te[im] + Te[i]); // electron temperature
    const BoutReal phi_wall = wall_potential[i];

    // Electron saturation at phi = phi_wall
    const BoutReal _potential = 0.5 * (phi[im] + phi[i]);
    const BoutReal phisheath = floor_potential ? floor(_potential, phi_wall) : _potential;

    // Electron velocity into sheath (< 0)
    // Equal to Bohm for single ions and no currents
    const BoutReal vesheath = bdry.toward * sqrt(tesheath / (TWOPI * Me)) * (1. - Ge)
                              * exp(-(phisheath - phi_wall) / floor(tesheath, 1e-5));

    const BoutReal vesheath_flow = no_flow ? 0.0 : vesheath;

    // Electron sheath heat transmission: either user-set or from model
    const BoutReal _gamma_e = [&]() {
      if constexpr (kind == hermes::SheathKind::simple) {
        return gamma_e;
      } else if constexpr (kind == hermes::SheathKind::normal) {
        return floor((2 / (1. - Ge)) + ((phisheath - phi_wall) / floor(tesheath, 1e-5)),
                     0.0);
      } else if constexpr (kind == hermes::SheathKind::insulating) {
        HERMES_UNREACHABLE;
      } else {
        static_assert(false, "Unhandled sheath kind");
      }
      HERMES_UNREACHABLE;
      return 0.0;
    }();

    const BoutReal adiabatic = -1 - (1 / (electron_adiabatic - 1));

    // Heat flux. Note: Here this is negative because vesheath < 0
    const BoutReal _q = (_gamma_e * tesheath * nesheath * vesheath)
                        // Take into account the flow of energy due to fluid flow
                        // This is additional energy flux through the sheath
                        + ((adiabatic * tesheath - 0.5 * Me * SQ(vesheath_flow))
                           * nesheath * vesheath_flow);

    const BoutReal q = bdry.is_lower() ? std::min(_q, 0.0) : std::max(_q, 0.0);

    Ve[im] = 2 * vesheath_flow - Ve[i];
    NVe[im] = 2 * Me * nesheath * vesheath_flow - NVe[i];

    // Cross-sectional area in XZ plane and cell volume
    const auto& da = (bdry.is_lower() ? da_m : da_p)[i];
    const auto& da_dv = (bdry.is_lower() ? da_dv_m : da_dv_p)[i];

    // Get power and energy source
    const BoutReal heatflow = q * da;        // [W]
    const BoutReal power = heatflow / dv[i]; // [Wm^-3]
    electron_energy_source[i] -= bdry.toward * power;

    // Total heat flux for diagnostic purposes
    const BoutReal q_no_flow = _gamma_e * tesheath * nesheath * vesheath_flow; // [Wm^-2]
    hflux_e[i] -= bdry.toward * q_no_flow * da_dv;                             // [Wm^-3]
    // lower Y: sheath boundary power placed in final domain cell
    // upper Y: sheath boundary power placed in ylow side of inner guard cell
    const auto i_power = bdry.is_lower() ? i : i.yp();
    electron_sheath_power_ylow[i_power] += heatflow; // [W]
  };

  if (lower_y) {
    loop_boundary<YBoundarySide::lower>(*mesh, set_electron_boundaries);
  }
  if (upper_y) {
    loop_boundary<YBoundarySide::upper>(*mesh, set_electron_boundaries);
  }

  // Set electron density and temperature, now with boundary conditions
  // Note: Clear parallel slices because they do not contain boundary conditions.
  Ne.clearParallelSlices();
  Te.clearParallelSlices();
  Pe.clearParallelSlices();
  setBoundary(electrons["density"], fromFieldAligned(Ne));
  setBoundary(electrons["temperature"], fromFieldAligned(Te));
  setBoundary(electrons["pressure"], fromFieldAligned(Pe));

  // Insulating BCs set these _after_ calculating ion flows
  if (not is_insulating()) {
    set(diagnostics["e"]["energy_source"], hflux_e);

    // Set energy source (negative in cell next to sheath)
    // Note: electron_energy_source includes any sources previously set in other
    // components
    set(electrons["energy_source"], fromFieldAligned(electron_energy_source));

    // Add the total sheath power flux to the tracker of y power flows
    add(electrons["energy_flow_ylow"], fromFieldAligned(electron_sheath_power_ylow));

    if (IS_SET_NOBOUNDARY(electrons["velocity"])) {
      Ve.clearParallelSlices();
      setBoundary(electrons["velocity"], fromFieldAligned(Ve));
    }
    if (IS_SET_NOBOUNDARY(electrons["momentum"])) {
      NVe.clearParallelSlices();
      setBoundary(electrons["momentum"], fromFieldAligned(NVe));
    }
  }

  if (always_set_phi or (state.isSection("fields") and state["fields"].isSet("phi"))) {
    // Set the potential, including boundary conditions
    phi.clearParallelSlices();
    setBoundary(state["fields"]["phi"], fromFieldAligned(phi));
  }

  //////////////////////////////////////////////////////////////////
  // Iterate through all ions
  for (const auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }

    Options& species = allspecies[kv.first]; // Note: Need non-const

    // Ion charge
    const BoutReal Zi = species.isSet("charge") ? get<BoutReal>(species["charge"]) : 0.0;

    if (Zi == 0.0) {
      continue; // Neutral -> skip
    }

    // Characteristics of this species
    const BoutReal Mi = get<BoutReal>(species["AA"]);
    // Ratio of specific heats (ideal gas)
    const BoutReal adiabatic =
        IS_SET(species["adiabatic"]) ? get<BoutReal>(species["adiabatic"]) : 5. / 3;

    // Density and temperature boundary conditions will be imposed (free)
    Field3D Ni = toFieldAligned(floor(getNoBoundary<Field3D>(species["density"]), 0.0));
    Field3D Ti = toFieldAligned(getNoBoundary<Field3D>(species["temperature"]));
    Field3D Pi = species.isSet("pressure")
                     ? toFieldAligned(getNoBoundary<Field3D>(species["pressure"]))
                     : Ni * Ti;

    // Get the velocity and momentum
    // These will be modified at the boundaries
    // and then put back into the state
    Field3D Vi = species.isSet("velocity")
                     ? toFieldAligned(getNoBoundary<Field3D>(species["velocity"]))
                     : zeroFrom(Ni);
    Field3D NVi = species.isSet("momentum")
                      ? toFieldAligned(getNoBoundary<Field3D>(species["momentum"]))
                      : Mi * Ni * Vi;

    // Energy source will be modified in the domain
    Field3D energy_source =
        species.isSet("energy_source")
            ? toFieldAligned(getNonFinal<Field3D>(species["energy_source"]))
            : zeroFrom(Ni);

    // Initialise sheath ion heat flux. This will be created for each species
    // saved in diagnostics struct and then destroyed and re-created for next species
    Field3D hflux_i = zeroFrom(energy_source);
    Field3D particle_source = zeroFrom(energy_source);
    // Field to capture total sheath heat flux for diagnostics
    Field3D ion_sheath_power_ylow = zeroFrom(Ne);

    const auto set_ion_boundary = [&](const YBoundary& bdry) {
      const auto i = bdry.i;
      const auto ip = bdry.ip;
      const auto im = bdry.im;

      // Free gradient of log electron density and temperature
      // This ensures that the guard cell values remain positive
      Ni[im] = limitFree(Ni[ip], Ni[i], density_boundary_mode);
      Ti[im] = limitFree(Ti[ip], Ti[i], temperature_boundary_mode);
      Pi[im] = limitFree(Pi[ip], Pi[i], pressure_boundary_mode);

      // Calculate sheath values at half-way points (cell edge)
      const BoutReal nesheath = 0.5 * (Ne[im] + Ne[i]);
      const BoutReal nisheath = 0.5 * (Ni[im] + Ni[i]);
      // electron temperature
      const BoutReal tesheath = floor(0.5 * (Te[im] + Te[i]), 1e-5);
      // ion temperature
      const BoutReal tisheath = floor(0.5 * (Ti[im] + Ti[i]), 1e-5);

      // Ion sheath heat transmission coefficient
      // Equation (22) in Tskhakaya 2005
      // with
      //
      // 1 / (1 + ∂_{ln n_e} ln s_i = s_i ∂_z n_e / ∂_z n_i
      // (from comparing C_i^2 in eq. 9 with eq. 20
      // Concentration
      const BoutReal s_i = std::clamp(nisheath / floor(nesheath, 1e-10), 0., 1.);
      BoutReal grad_ne = Ne[i] - nesheath;
      BoutReal grad_ni = Ni[i] - nisheath;

      if (std::abs(grad_ni) < 1.e-3) {
        grad_ni = grad_ne = 1e-3; // Remove kinetic correction term
      }

      // Ion speed into sheath
      const BoutReal C_i_sq = [&]() {
        if constexpr (kind == hermes::SheathKind::normal
                      or kind == hermes::SheathKind::insulating) {
          return ((adiabatic * tisheath) + (Zi * s_i * tesheath) * (grad_ne / grad_ni))
                 / Mi;
        } else if constexpr (kind == hermes::SheathKind::simple) {
          return ((sheath_ion_polytropic * tisheath) + (Zi * tesheath)) / Mi;
        } else {
          static_assert(false, "Unhandled sheath kind");
        }
      }();

      // Negative -> into sheath
      const BoutReal visheath = bdry.is_lower() ? std::min(Vi[i], -sqrt(C_i_sq))
                                                : std::max(Vi[i], sqrt(C_i_sq));
      const BoutReal visheath_flow = no_flow ? 0.0 : visheath;

      // Note: Here this is negative because visheath < 0
      const BoutReal q = [&]() {
        if constexpr (kind == hermes::SheathKind::normal
                      or kind == hermes::SheathKind::insulating) {
          // Ion sheath heat transmission coefficient
          // TODO(peter): Should this 2.5 be `-1 - 1/(adiabatic - 1)`?
          const BoutReal _gamma_i = 2.5 + (0.5 * Mi * C_i_sq / tisheath);

          const auto _q =
              ((_gamma_i - 1 - 1 / (adiabatic - 1)) * tisheath - 0.5 * Mi * C_i_sq)
              * nisheath * visheath;
          return bdry.is_lower() ? std::min(_q, 0.0) : std::max(_q, 0.0);
        } else if constexpr (kind == hermes::SheathKind::simple) {
          return (gamma_i * tisheath * nisheath * visheath)
                 // Take into account the flow of energy due to fluid flow
                 // This is additional energy flux through the sheath
                 - ((2.5 * tisheath + 0.5 * Mi * SQ(visheath_flow)) * nisheath
                    * visheath_flow);
        } else {
          static_assert(false, "Unhandled sheath kind");
        }
      }();

      // Set boundary conditions on flows
      Vi[im] = 2. * visheath_flow - Vi[i];
      NVi[im] = 2. * Mi * nisheath * visheath_flow - NVi[i];

      if constexpr (is_insulating()) {
        // Add electron flow to balance current
        Ve[im] += 2. * visheath * Zi;
        NVe[im] += 2. * Me * nisheath * visheath;
      }

      // Cross-sectional area in XZ plane and cell volume
      const auto& da = (bdry.is_lower() ? da_m : da_p)[i];
      const auto& da_dv = (bdry.is_lower() ? da_dv_m : da_dv_p)[i];

      // Get power and energy source
      // Multiply by cell area to get power
      const BoutReal heatflow = q * da; // [W]
      // Divide by cell volume to get energy loss rate
      const BoutReal power = heatflow / dv[i]; // [Wm^-3]
      ASSERT2(std::isfinite(power));
      energy_source[i] -= bdry.toward * power;
      // Diagnostics only
      particle_source[i] -= nisheath * visheath_flow * da_dv; // [m^-3s^-1]

      // Total heat flux for diagnostic purposes
      const BoutReal q_no_flow = gamma_i * tisheath * nisheath * visheath_flow; // [Wm^-2]
      hflux_i[i] -= bdry.toward * q_no_flow * da_dv;                            // [Wm^-3]
      // lower Y: sheath boundary power placed in final domain cell
      // upper Y: sheath boundary power placed in ylow side of inner guard cell
      const auto i_power = bdry.is_lower() ? i : i.yp();
      ion_sheath_power_ylow[i_power] += heatflow;
    };

    if (lower_y) {
      loop_boundary<YBoundarySide::lower>(*mesh, set_ion_boundary);
    }
    if (upper_y) {
      loop_boundary<YBoundarySide::upper>(*mesh, set_ion_boundary);
    }

    // Finished boundary conditions for this species
    // Put the modified fields back into the state.

    Ni.clearParallelSlices();
    Ti.clearParallelSlices();
    Pi.clearParallelSlices();
    setBoundary(species["density"], fromFieldAligned(Ni));
    setBoundary(species["temperature"], fromFieldAligned(Ti));
    setBoundary(species["pressure"], fromFieldAligned(Pi));

    if (species.isSet("velocity")) {
      Vi.clearParallelSlices();
      setBoundary(species["velocity"], fromFieldAligned(Vi));
    }

    if (species.isSet("momentum")) {
      NVi.clearParallelSlices();
      setBoundary(species["momentum"], fromFieldAligned(NVi));
    }

    // Additional loss of energy through sheath
    // Note: energy_source already includes previously set values
    set(species["energy_source"], fromFieldAligned(energy_source));

    // Add the total sheath power flux to the tracker of y power flows
    add(species["energy_flow_ylow"], fromFieldAligned(ion_sheath_power_ylow));

    set(diagnostics[species.name()]["energy_source"], hflux_i);
    set(diagnostics[species.name()]["particle_source"], particle_source);
  }

  if (is_insulating()) {
    //////////////////////////////////////////////////////////////////
    // Electrons
    // This time adding energy sink, having calculated flow

    const auto set_electron_source = [&](const YBoundary& bdry) {
      const auto i = bdry.i;
      const auto ip = bdry.ip;
      const auto im = bdry.im;

      const BoutReal nesheath = 0.5 * (Ne[im] + Ne[i]);
      const BoutReal tesheath = 0.5 * (Te[im] + Te[i]); // electron temperature
      // Electron velocity into sheath (< 0). Calculated from ion flow
      const BoutReal vesheath = 0.5 * (Ve[im] + Ve[i]);
      const BoutReal vesheath_flow = no_flow ? 0.0 : vesheath;

      // Take into account the flow of energy due to fluid flow
      // This is additional energy flux through the sheath
      // Note: Here this is negative because vesheath < 0
      BoutReal q = ((gamma_e - 1 - 1 / (electron_adiabatic - 1)) * tesheath
                    - 0.5 * Me * SQ(vesheath))
                   * nesheath * vesheath;

      // Cross-sectional area in XZ plane and cell volume
      const auto& da = (bdry.is_lower() ? da_m : da_p)[i];
      const auto& da_dv = (bdry.is_lower() ? da_dv_m : da_dv_p)[i];

      // Get power and energy source
      // Multiply by cell area to get power
      const BoutReal heatflow = q * da; // [W]
      // Divide by cell volume to get energy loss rate
      const BoutReal power = heatflow / dv[i]; // [Wm^-3]

#if CHECKLEVEL >= 1
      if (!std::isfinite(power)) {
        throw BoutException("Non-finite power at {} : Te {} Ne {} Ve {}", i, tesheath,
                            nesheath, vesheath);
      }
#endif

      electron_energy_source[i] -= bdry.toward * power;
      // Total heat flux for diagnostic purposes
      const BoutReal q_no_flow = gamma_e * tesheath * nesheath * vesheath_flow; // [Wm^-2]
      hflux_e[i] -= bdry.toward * q_no_flow * da_dv;                            // [Wm^-3]
      // lower Y: sheath boundary power placed in final domain cell
      // upper Y: sheath boundary power placed in ylow side of inner guard cell
      const auto i_power = bdry.is_lower() ? i : i.yp();
      electron_sheath_power_ylow[i_power] += heatflow; // [W]
    };

    if (lower_y) {
      loop_boundary<YBoundarySide::lower>(*Ne.getMesh(), set_electron_source);
    }
    if (upper_y) {
      loop_boundary<YBoundarySide::upper>(*Ne.getMesh(), set_electron_source);
    }

    set(diagnostics["e"]["energy_source"], hflux_e);

    // Set energy source (negative in cell next to sheath)
    // Note: electron_energy_source includes any sources previously set in other
    // components
    set(electrons["energy_source"], fromFieldAligned(electron_energy_source));

    // Add the total sheath power flux to the tracker of y power flows
    add(electrons["energy_flow_ylow"], fromFieldAligned(electron_sheath_power_ylow));

    if (IS_SET_NOBOUNDARY(electrons["velocity"])) {
      Ve.clearParallelSlices();
      setBoundary(electrons["velocity"], fromFieldAligned(Ve));
    }
    if (IS_SET_NOBOUNDARY(electrons["momentum"])) {
      NVe.clearParallelSlices();
      setBoundary(electrons["momentum"], fromFieldAligned(NVe));
    }
  }
}

template <hermes::SheathKind kind>
void SheathBoundaryBase<kind>::outputVars(Options& state) {
  AUTO_TRACE();
  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

  if (diagnose) {
    /// Iterate through the first species in each collision pair
    const std::map<std::string, Options>& level1 = diagnostics.getChildren();
    for (const auto& [species_name, section] : level1) {
      set_with_attrs(state[fmt::format("E{}_sheath", species_name)],
                     getNonFinal<Field3D>(section["energy_source"]),
                     {{"time_dimension", "t"},
                      {"units", "W / m^3"},
                      {"conversion", Pnorm * Omega_ci},
                      {"standard_name", "energy source"},
                      {"long_name", species_name + " sheath energy source"},
                      {"source", "sheath_boundary_simple"}});

      if (species_name != "e") {
        set_with_attrs(state[fmt::format("S{}_sheath", species_name)],
                       getNonFinal<Field3D>(section["particle_source"]),
                       {{"time_dimension", "t"},
                        {"units", "m^-3 s^-1"},
                        {"conversion", Nnorm * Omega_ci},
                        {"standard_name", "energy source"},
                        {"long_name", species_name + " sheath energy source"},
                        {"source", "sheath_boundary_simple"}});
      }
    }
  }
}
