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
BoutReal clip(BoutReal value, BoutReal min, BoutReal max) {
  if (value < min)
    return min;
  if (value > max)
    return max;
  return value;
}


// stubs to turn to/from field aligned into no-ops for FCI
template <typename F>
F fromFieldAlignedAuto(const F& f) {
  if (f.isFci()) {
    return f;
  }
  return ::fromFieldAligned(f);
}
template <typename F>
F toFieldAlignedAuto(const F& f) {
  if (f.isFci()) {
    return f;
  }
  return ::toFieldAligned(f);
}

/// Calculate potential phi assuming zero current
Field3D calculate_phi_simple(Options& allspecies, YBoundary yboundary, const Field3D& Ne,
                             const Field3D& Te, BoutReal Me, BoutReal Ge,
                             const Field3D& wall_potential,
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
                     ? toFieldAlignedAuto(getNoBoundary<Field3D>(species["velocity"]))
                     : zeroFrom(Ni);

    yboundary.iter_pnts([&](auto& pnt) {
      // Free gradient of log density and temperature
      // This ensures that the guard cell values remain positive

      // Calculate sheath values at half-way points (cell edge)
      const BoutReal nisheath = pnt.extrapolate_sheath_free(Ni, density_boundary_mode);
      // electron temperature
      const BoutReal tesheath = floor(pnt.extrapolate_sheath_free(Te, temperature_boundary_mode), 1e-5);
      // ion temperature
      const BoutReal tisheath = floor(pnt.extrapolate_sheath_free(Ti, temperature_boundary_mode), 1e-5);

      // Sound speed
      const BoutReal C_i = sqrt((sheath_ion_polytropic * tisheath + Zi * tesheath) / Mi);
      const BoutReal visheath =
          pnt.dir < 0 ? std::min(pnt.ythis(Vi), -C_i) : std::max(pnt.ythis(Vi), C_i);

      if (pnt.abs_offset() == 1) {
        pnt.ythis(ion_sum) += pnt.dir * Zi * nisheath * visheath;
      }
    });
  }

  Field3D phi = emptyFrom(Ne);

  yboundary.iter_pnts([&](auto& pnt) {
    if (pnt.abs_offset() == 1) {
      auto i = pnt.ind();
      // Calculate sheath values at cell edge
      const BoutReal nesheath = pnt.extrapolate_sheath_free(Ne, density_boundary_mode);
      const BoutReal tesheath = floor(pnt.extrapolate_sheath_free(Te, temperature_boundary_mode), 1e-5);
      
      phi[i] =
        tesheath * log(sqrt(tesheath / (Me * TWOPI)) * (1. - Ge) * nesheath / ion_sum[i]);

      // Add bias potential
      phi[i] += wall_potential[i];
      // Constant into sheath
      pnt.neumann_o1(phi, 0);
    }
  });

  return phi;
}

/// Calculate potential phi assuming zero current.
/// Note: This is equation (22) in Tskhakaya 2005, with I = 0
Field3D calculate_phi_Tskhakaya(Options& allspecies, YBoundary yboundary,
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

    const Field3D Ni = toFieldAlignedAuto(
        floor(GET_NOBOUNDARY(Field3D, species["density"]).asField3DParallel(), 0.0));
    const Field3D Ti =
        toFieldAlignedAuto(GET_NOBOUNDARY(Field3D, species["temperature"]));
    const BoutReal Mi = GET_NOBOUNDARY(BoutReal, species["AA"]);
    const BoutReal Zi = GET_NOBOUNDARY(BoutReal, species["charge"]);

    const BoutReal adiabatic = IS_SET(species["adiabatic"])
                                   ? get<BoutReal>(species["adiabatic"])
                                   : 5. / 3; // Ratio of specific heats (ideal gas)

    yboundary.iter_pnts([&](auto& pnt) {
      if (pnt.abs_offset() == 1) {
	// Free boundary extrapolate ion concentration
	// Limit range to [0,1]
	BoutReal s_i =
	  clip(pnt.extrapolate_sheath_o2([&, Ni, Ne](int yoffset, Ind3D ind) {
	    return Ni.ynext(yoffset)[ind] / Ne.ynext(yoffset)[ind];
	  }),
	    0.0, 1.0);

	if (!std::isfinite(s_i)) {
	  s_i = 1.0;
	}
	const BoutReal te = pnt.ythis(Te);
	const BoutReal ti = pnt.ythis(Ti);

	// Equation (9) in Tskhakaya 2005
	BoutReal grad_ne = -pnt.extrapolate_grad_o2(Ne);
	BoutReal grad_ni = -pnt.extrapolate_grad_o2(Ni);

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
        ion_sum[pnt.ind()] += s_i * Zi * sin_alpha * sqrt(C_i_sq);
      }
    });
  }

  Field3D phi = emptyFrom(Ne.asField3DParallel()); // So phi is field aligned

  yboundary.iter_pnts([&](auto& pnt) {
    if (pnt.abs_offset() == 1) {
      const auto i = pnt.ind();
      if (Te[i] <= 0.0) {
	phi[i] = 0.0;
      } else {
	phi[i] = Te[i] * log(sqrt(Te[i] / (Me * TWOPI)) * (1. - Ge) / ion_sum[i]);
      }
      phi[i] += wall_potential[i];        // Add bias potential
      pnt.setAll(phi, phi[i]);            // Constant into sheath
    }
  });

  return phi;
}

// This is a workaround before CWG2518/P2593R1, taken from cppreference.com
template<hermes::SheathKind>
constexpr bool dependent_false = false;

/// If phi is set, use free boundary condition;
/// If phi not set, calculate assuming zero current, according to sheath kind
template <hermes::SheathKind kind>
Field3D
init_phi(Options& state, Options& allspecies, YBoundary yboundary, const Field3D& Ne,
         const Field3D& Te, BoutReal Me, BoutReal Ge, const Field3D& wall_potential,
         BoutReal sin_alpha, SheathLimitMode density_boundary_mode,
         SheathLimitMode temperature_boundary_mode, BoutReal sheath_ion_polytropic) {

  if (IS_SET_NOBOUNDARY(state["fields"]["phi"])) {
    return toFieldAlignedAuto(getNoBoundary<Field3D>(state["fields"]["phi"]));
  }

  if constexpr (kind == hermes::SheathKind::normal) {
    return calculate_phi_Tskhakaya(allspecies, yboundary, Ne, Te, Me, Ge, wall_potential,
                                   sin_alpha);
  } else if constexpr (kind == hermes::SheathKind::simple) {
    return calculate_phi_simple(allspecies, yboundary, Ne, Te, Me, Ge, wall_potential,
                                density_boundary_mode, temperature_boundary_mode,
                                sheath_ion_polytropic);
  } else if constexpr (kind == hermes::SheathKind::insulating) {
    // No model for potential for this case, but return something to
    // avoid duplicated checks later
    return 0.0;
  } else {
    static_assert(dependent_false<kind>, "Unhandled sheath kind");
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

  // init parallel bc iterator
  yboundary.init(options);

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
  if (wall_potential.isFci()) {
    Coordinates* coord = mesh->getCoordinates();
    if (!coord->dx.hasParallelSlices()) {
      mesh->communicate(coord->dx);
      coord->dx.applyParallelBoundary("parallel_neumann_o2");
    }
    if (!coord->dz.hasParallelSlices()) {
      mesh->communicate(coord->dz);
      coord->dz.applyParallelBoundary("parallel_neumann_o2");
    }
  }
}

template <hermes::SheathKind kind>
void SheathBoundaryBase<kind>::transform(Options& state) {
  AUTO_TRACE();

  Options& allspecies = state["species"];
  Options& electrons = allspecies["e"];

  // Need electron properties
  // Not const because boundary conditions will be set
  Field3D Ne = toFieldAligned(
      floor(GET_NOBOUNDARY(Field3D, electrons["density"]).asField3DParallel(), 0.0));
  Field3D Te = toFieldAligned(GET_NOBOUNDARY(Field3D, electrons["temperature"]));
  Field3D Pe = IS_SET_NOBOUNDARY(electrons["pressure"])
                   ? toFieldAligned(getNoBoundary<Field3D>(electrons["pressure"]))
                   : Te * Ne.asField3DParallel();

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
                   : zeroFrom(Ne.asField3DParallel());

  Field3D NVe = IS_SET_NOBOUNDARY(electrons["momentum"])
                    ? toFieldAligned(getNoBoundary<Field3D>(electrons["momentum"]))
                    : zeroFrom(Ne.asField3DParallel());

  Coordinates* coord = mesh->getCoordinates();

  //////////////////////////////////////////////////////////////////
  // Electrostatic potential
  Field3D phi = init_phi<kind>(state, allspecies, yboundary, Ne, Te, Me, Ge,
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

  yboundary.iter_pnts([&](auto& pnt) {
    if (pnt.abs_offset() != 1) {
      return;
    }
    auto i = pnt.ind();
    auto da_c = pnt.dir == 1 ? da_p : da_m;
    da_c[i] = pnt.interpolate_sheath_o2(coord->J)
              / pnt.interpolate_sheath_o2([&, coord](int yoffset, Ind3D i) {
                  return sqrt(coord->g_22.ynext(yoffset)[i]);
                })
              * pnt.interpolate_sheath_o2(coord->dx)
              * pnt.interpolate_sheath_o2(coord->dz); // [m^2]

    dv[i] = (coord->dx[i] * coord->dy[i] * coord->dz[i] * coord->J[i]);

    (pnt.dir == 1 ? da_dv_p : da_dv_m)[i] = da_c[i] / dv[i];
  });

  yboundary.iter_pnts([&](auto& pnt) {
    // Free gradient of log electron density and temperature
    // Limited so that the values don't increase into the sheath
    // This ensures that the guard cell values remain positive
    pnt.set_free(Ne, density_boundary_mode);
    pnt.set_free(Te, temperature_boundary_mode);
    pnt.set_free(Pe, pressure_boundary_mode);

    // Free boundary potential linearly extrapolated
    pnt.set_free(phi, SheathLimitMode::linear_free);

    if constexpr (is_insulating()) {
      // Set zero flow through boundary
      // This will be modified when iterating over the ions
      pnt.dirichlet_o2(Ve, 0);
      pnt.dirichlet_o2(NVe, 0);
      return;
    }

    const BoutReal nesheath = pnt.interpolate_sheath_o2(Ne);
    const BoutReal tesheath = pnt.interpolate_sheath_o2(Te); // electron temperature
    const BoutReal phi_wall = wall_potential[pnt.ind()];

    // Electron saturation at phi = phi_wall
    const BoutReal _potential = pnt.interpolate_sheath_o2(phi);
    const BoutReal phisheath = floor_potential ? floor(_potential, phi_wall) : _potential;

    // Electron velocity into sheath (< 0)
    // Equal to Bohm for single ions and no currents
    const BoutReal vesheath = pnt.dir * sqrt(tesheath / (TWOPI * Me)) * (1. - Ge)
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
        static_assert(dependent_false<kind>, "Unhandled sheath kind");
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

    const BoutReal q = pnt.dir == -1 ? std::min(_q, 0.0) : std::max(_q, 0.0);

    pnt.dirichlet_o2(Ve, vesheath_flow);
    pnt.dirichlet_o2(NVe, Me * nesheath * vesheath_flow);

    if (pnt.abs_offset() == 1) {
      auto i = pnt.ind();
      // Cross-sectional area in XZ plane and cell volume
      const auto& da = (pnt.is_lower() ? da_m : da_p)[i];
      const auto& da_dv = (pnt.is_lower() ? da_dv_m : da_dv_p)[i];

      // Get power and energy source
      const BoutReal heatflow = q * da;        // [W]
      const BoutReal power = heatflow / dv[i]; // [Wm^-3]
      electron_energy_source[i] -= pnt.dir * power;

      // Total heat flux for diagnostic purposes
      const BoutReal q_no_flow = _gamma_e * tesheath * nesheath * vesheath_flow; // [Wm^-2]
      hflux_e[i] -= pnt.dir * q_no_flow * da_dv; // [Wm^-3]
      // TODO: Fix for FCI
      // lower Y: sheath boundary power placed in final domain cell
      // upper Y: sheath boundary power placed in ylow side of inner guard cell
      const auto i_power = pnt.is_lower() ? i : i.yp();
      electron_sheath_power_ylow[i_power] += heatflow; // [W]
    }
  });

  // Set electron density and temperature, now with boundary conditions
  // Note: Clear parallel slices because they do not contain boundary conditions.
  if (! Ne.isFci()) {
    Ne.clearParallelSlices();
    Te.clearParallelSlices();
    Pe.clearParallelSlices();
  }
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
      if (! Ve.isFci()) {
	Ve.clearParallelSlices();
      }
      setBoundary(electrons["velocity"], fromFieldAligned(Ve));
    }
    if (IS_SET_NOBOUNDARY(electrons["momentum"])) {
      if (! NVe.isFci()) {
	NVe.clearParallelSlices();
      }
      setBoundary(electrons["momentum"], fromFieldAligned(NVe));
    }
  }

  if (always_set_phi or (state.isSection("fields") and state["fields"].isSet("phi"))) {
    // Set the potential, including boundary conditions
    if (! phi.isFci()) {
      phi.clearParallelSlices();
    }
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
    Field3D Ni = toFieldAligned(floor(getNoBoundary<Field3D>(species["density"]).asField3DParallel(), 0.0));
    Field3D Ti = toFieldAligned(getNoBoundary<Field3D>(species["temperature"]));
    Field3D Pi = species.isSet("pressure")
                     ? toFieldAligned(getNoBoundary<Field3D>(species["pressure"]))
      : Ni * Ti.asField3DParallel();

    // Get the velocity and momentum
    // These will be modified at the boundaries
    // and then put back into the state
    Field3D Vi = species.isSet("velocity")
                     ? toFieldAligned(getNoBoundary<Field3D>(species["velocity"]))
      : zeroFrom(Ni.asField3DParallel());
    Field3D NVi = species.isSet("momentum")
                      ? toFieldAligned(getNoBoundary<Field3D>(species["momentum"]))
      : Mi * Ni * Vi.asField3DParallel();

    // Energy source will be modified in the domain
    Field3D energy_source =
        species.isSet("energy_source")
            ? toFieldAligned(getNonFinal<Field3D>(species["energy_source"]))
            : zeroFrom(Ni.asField3DParallel());

    // Initialise sheath ion heat flux. This will be created for each species
    // saved in diagnostics struct and then destroyed and re-created for next species
    Field3D hflux_i = zeroFrom(energy_source);
    Field3D particle_source = zeroFrom(energy_source);
    // Field to capture total sheath heat flux for diagnostics
    Field3D ion_sheath_power_ylow = zeroFrom(Ne);

    yboundary.iter_pnts([&](auto& pnt) {
      // Free gradient of log electron density and temperature
      // This ensures that the guard cell values remain positive
      pnt.set_free(Ni, density_boundary_mode);
      pnt.set_free(Ti, temperature_boundary_mode);
      pnt.set_free(Pi, pressure_boundary_mode);

      // Calculate sheath values at half-way points (cell edge)
      const BoutReal nesheath = pnt.interpolate_sheath_o2(Ne);
      const BoutReal nisheath = pnt.interpolate_sheath_o2(Ni);
      // electron temperature
      const BoutReal tesheath = floor(pnt.interpolate_sheath_o2(Te), 1e-5);
      // ion temperature
      const BoutReal tisheath = floor(pnt.interpolate_sheath_o2(Ti), 1e-5);

      // Ion sheath heat transmission coefficient
      // Equation (22) in Tskhakaya 2005
      // with
      //
      // 1 / (1 + ∂_{ln n_e} ln s_i = s_i ∂_z n_e / ∂_z n_i
      // (from comparing C_i^2 in eq. 9 with eq. 20
      // Concentration
      const BoutReal s_i = std::clamp(nisheath / floor(nesheath, 1e-10), 0., 1.);
      BoutReal grad_ne = pnt.extrapolate_grad_o2(Ne) * 2;
      BoutReal grad_ni = pnt.extrapolate_grad_o2(Ni) * 2;

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
          static_assert(dependent_false<kind>, "Unhandled sheath kind");
        }
      }();

      // Negative -> into sheath
      const BoutReal visheath = pnt.is_lower() ? std::min(pnt.ythis(Vi), -sqrt(C_i_sq))
                                               : std::max(pnt.ythis(Vi), sqrt(C_i_sq));
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
          return pnt.is_lower() ? std::min(_q, 0.0) : std::max(_q, 0.0);
        } else if constexpr (kind == hermes::SheathKind::simple) {
          return (gamma_i * tisheath * nisheath * visheath)
                 // Take into account the flow of energy due to fluid flow
                 // This is additional energy flux through the sheath
                 - ((2.5 * tisheath + 0.5 * Mi * SQ(visheath_flow)) * nisheath
                    * visheath_flow);
        } else {
          static_assert(dependent_false<kind>, "Unhandled sheath kind");
        }
      }();

      // Set boundary conditions on flows
      pnt.dirichlet_o2(Vi, visheath_flow);
      pnt.dirichlet_o2(NVi, Mi * nisheath * visheath_flow);

      if constexpr (is_insulating()) {
        // Add electron flow to balance current
	if (Ve.isFci()) {
	  throw BoutException("insulation not yet implemented for FCI");
        }
        auto im = pnt.ind().yp(pnt.dir);
        Ve[im] += 2. * visheath * Zi;
        NVe[im] += 2. * Me * nisheath * visheath;
      }

      if (pnt.abs_offset() == 1) {
        auto i = pnt.ind();
        // Cross-sectional area in XZ plane and cell volume
        const auto& da = (pnt.is_lower() ? da_m : da_p)[i];
        const auto& da_dv = (pnt.is_lower() ? da_dv_m : da_dv_p)[i];

        // Get power and energy source
        // Multiply by cell area to get power
        const BoutReal heatflow = q * da; // [W]
	// Divide by cell volume to get energy loss rate
	const BoutReal power = heatflow / dv[i]; // [Wm^-3]
	ASSERT2(std::isfinite(power));
        energy_source[i] -= pnt.dir * power;
        // Diagnostics only
        particle_source[i] -= nisheath * visheath_flow * da_dv; // [m^-3s^-1]

	// Total heat flux for diagnostic purposes
	const BoutReal q_no_flow = gamma_i * tisheath * nisheath * visheath_flow; // [Wm^-2]
        hflux_i[i] -= pnt.dir * q_no_flow * da_dv; // [Wm^-3]
        // TODO: Fix for FCI
        // lower Y: sheath boundary power placed in final domain cell
        // upper Y: sheath boundary power placed in ylow side of inner guard cell
        const auto i_power = pnt.is_lower() ? i : i.yp();
        ion_sheath_power_ylow[i_power] += heatflow;
      }
    });

    // Finished boundary conditions for this species
    // Put the modified fields back into the state.

    if (! Ni.isFci()) {
      Ni.clearParallelSlices();
      Ti.clearParallelSlices();
      Pi.clearParallelSlices();
    }
    setBoundary(species["density"], fromFieldAligned(Ni));
    setBoundary(species["temperature"], fromFieldAligned(Ti));
    setBoundary(species["pressure"], fromFieldAligned(Pi));

    if (species.isSet("velocity")) {
      if (! Vi.isFci()) {
	Vi.clearParallelSlices();
      }
      setBoundary(species["velocity"], fromFieldAligned(Vi));
    }

    if (species.isSet("momentum")) {
      if (! NVi.isFci()) {
	NVi.clearParallelSlices();
      }
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

    yboundary.iter_pnts([&](auto& pnt) {
      const BoutReal nesheath = pnt.interpolate_sheath_o2(Ne);
      const BoutReal tesheath = pnt.interpolate_sheath_o2(Te);

      // Electron velocity into sheath (< 0). Calculated from ion flow
      const BoutReal vesheath = pnt.interpolate_sheath_o2(Ve);
      const BoutReal vesheath_flow = no_flow ? 0.0 : vesheath;

      // Take into account the flow of energy due to fluid flow
      // This is additional energy flux through the sheath
      // Note: Here this is negative because vesheath < 0
      BoutReal q = ((gamma_e - 1 - 1 / (electron_adiabatic - 1)) * tesheath
                    - 0.5 * Me * SQ(vesheath))
                   * nesheath * vesheath;

      if (pnt.abs_offset() == 1) {
	const auto i = pnt.ind();
	// Cross-sectional area in XZ plane and cell volume
        const auto& da = (pnt.is_lower() ? da_m : da_p)[i];
        const auto& da_dv = (pnt.is_lower() ? da_dv_m : da_dv_p)[i];

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

        electron_energy_source[i] -= pnt.dir * power;
        // Total heat flux for diagnostic purposes
        const BoutReal q_no_flow = gamma_e * tesheath * nesheath * vesheath_flow; // [Wm^-2]
        hflux_e[i] -= pnt.dir * q_no_flow * da_dv; // [Wm^-3]
        // lower Y: sheath boundary power placed in final domain cell
        // upper Y: sheath boundary power placed in ylow side of inner guard cell
        const auto i_power = pnt.is_lower() ? i : i.yp();
        electron_sheath_power_ylow[i_power] += heatflow; // [W]
      }
    });

    set(diagnostics["e"]["energy_source"], hflux_e);

    // Set energy source (negative in cell next to sheath)
    // Note: electron_energy_source includes any sources previously set in other
    // components
    set(electrons["energy_source"], fromFieldAligned(electron_energy_source));

    // Add the total sheath power flux to the tracker of y power flows
    add(electrons["energy_flow_ylow"], fromFieldAligned(electron_sheath_power_ylow));

    if (IS_SET_NOBOUNDARY(electrons["velocity"])) {
      if (! Ve.isFci()) {
	Ve.clearParallelSlices();
      }
      setBoundary(electrons["velocity"], fromFieldAligned(Ve));
    }
    if (IS_SET_NOBOUNDARY(electrons["momentum"])) {
      if (! NVe.isFci()) {
	NVe.clearParallelSlices();
      }
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
