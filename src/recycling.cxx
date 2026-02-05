
#include "../include/recycling.hxx"
#include "../include/hermes_utils.hxx" // For indexAt

#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/mesh.hxx>
#include <bout/utils.hxx> // for trim, strsplit

#include <algorithm>

using bout::globals::mesh;

Recycling::Recycling(std::string name, Options& alloptions, Solver*)
    : Component({readOnly("species:{from}:{from_inputs}"),
                 readOnly("species:{to}:{to_inputs}"),
                 readWrite("species:{to}:{outputs}")}) {

  const Options& units = alloptions["units"];
  const BoutReal Tnorm = units["eV"];

  Options& options = alloptions[name];

  auto species_list = strsplit(options["species"]
                                   .doc("Comma-separated list of species to recycle")
                                   .as<std::string>(),
                               ',');

  std::set<std::string> from_species, to_species;

  // Neutral pump
  // Mark cells as having a pump by setting the Field2D is_pump to 1 in the grid file
  // Works only on SOL and PFR edges, where it locally modifies the recycle multiplier to
  // the pump albedo
  is_pump = 0.0;
  mesh->get(is_pump, std::string("is_pump"));

  // We need sheath heat flux for fast recycling to work. We can't get it because
  // the advection component is in the finally() of evolve_pressure, so we must
  // recalculate it. This is only possible for sheath_boundary_simple where
  // gamma is known.
  if (alloptions.isSection("sheath_boundary_simple")) {
    has_sheath_boundary_simple = true;
    if (alloptions["sheath_boundary_simple"]["gamma_i"].isSet()) {
      gamma_i = alloptions["sheath_boundary_simple"]["gamma_i"];
    } else {
      gamma_i = 3.5; // Default value
    }
  } else {
    has_sheath_boundary_simple = false;
    gamma_i = 0.0;
  }

  for (const auto& species : species_list) {
    std::string from = trim(species, " \t\r()"); // The species name in the list

    if (from.empty())
      continue; // Missing

    // Get the options for this species
    Options& from_options = alloptions[from];
    std::string to = from_options["recycle_as"]
                         .doc("Name of the species to recycle into")
                         .as<std::string>();

    from_species.insert(from);
    to_species.insert(to);

    density_floor =
        options["density_floor"].doc("Minimum density floor").withDefault(1e-7);
    pressure_floor = density_floor * (1. / get<BoutReal>(alloptions["units"]["eV"]));

    diagnose = from_options["diagnose"]
                   .doc("Save additional diagnostics?")
                   .withDefault<bool>(false);

    BoutReal target_recycle_multiplier = from_options["target_recycle_multiplier"]
                                             .doc("Multiply the target recycled flux by "
                                                  "this factor. Should be >=0 and <= 1")
                                             .withDefault<BoutReal>(1.0);

    BoutReal sol_recycle_multiplier =
        from_options["sol_recycle_multiplier"]
            .doc("Multiply the sol recycled flux by this factor. Should be >=0 and <= 1")
            .withDefault<BoutReal>(1.0);

    BoutReal pfr_recycle_multiplier =
        from_options["pfr_recycle_multiplier"]
            .doc("Multiply the pfr recycled flux by this factor. Should be >=0 and <= 1")
            .withDefault<BoutReal>(1.0);

    BoutReal pump_recycle_multiplier =
        from_options["pump_recycle_multiplier"]
            .doc("Multiply the pump boundary recycling flux by this factor (like "
                 "albedo). Should be >=0 and <= 1")
            .withDefault<BoutReal>(1.0);

    BoutReal target_recycle_energy =
        from_options["target_recycle_energy"]
            .doc("Fixed energy of the recycled particles at target [eV]")
            .withDefault<BoutReal>(3.0)
        / Tnorm; // Normalise from eV

    BoutReal sol_recycle_energy =
        from_options["sol_recycle_energy"]
            .doc("Fixed energy of the recycled particles at sol [eV]")
            .withDefault<BoutReal>(3.0)
        / Tnorm; // Normalise from eV

    BoutReal pfr_recycle_energy =
        from_options["pfr_recycle_energy"]
            .doc("Fixed energy of the recycled particles at pfr [eV]")
            .withDefault<BoutReal>(3.0)
        / Tnorm; // Normalise from eV

    BoutReal target_fast_recycle_fraction =
        from_options["target_fast_recycle_fraction"]
            .doc("Fraction of ions undergoing fast reflection at target")
            .withDefault<BoutReal>(0);

    BoutReal pfr_fast_recycle_fraction =
        from_options["pfr_fast_recycle_fraction"]
            .doc("Fraction of ions undergoing fast reflection at pfr")
            .withDefault<BoutReal>(0);

    BoutReal sol_fast_recycle_fraction =
        from_options["sol_fast_recycle_fraction"]
            .doc("Fraction of ions undergoing fast reflection at sol")
            .withDefault<BoutReal>(0);

    BoutReal target_fast_recycle_energy_factor =
        from_options["target_fast_recycle_energy_factor"]
            .doc("Fraction of energy retained by fast recycled neutrals at target")
            .withDefault<BoutReal>(0);

    BoutReal sol_fast_recycle_energy_factor =
        from_options["sol_fast_recycle_energy_factor"]
            .doc("Fraction of energy retained by fast recycled neutrals at sol")
            .withDefault<BoutReal>(0);

    BoutReal pfr_fast_recycle_energy_factor =
        from_options["pfr_fast_recycle_energy_factor"]
            .doc("Fraction of energy retained by fast recycled neutrals at pfr")
            .withDefault<BoutReal>(0);

    if ((target_recycle_multiplier < 0.0) or (target_recycle_multiplier > 1.0)
        or (sol_recycle_multiplier < 0.0) or (sol_recycle_multiplier > 1.0)
        or (pfr_recycle_multiplier < 0.0) or (pfr_recycle_multiplier > 1.0)
        or (pump_recycle_multiplier < 0.0) or (pump_recycle_multiplier > 1.0)) {
      throw BoutException("All recycle multipliers must be betweeen 0 and 1");
    }

    // Populate recycling channel vector
    channels.push_back({from, to, target_recycle_multiplier, sol_recycle_multiplier,
                        pfr_recycle_multiplier, pump_recycle_multiplier,
                        target_recycle_energy, sol_recycle_energy, pfr_recycle_energy,
                        target_fast_recycle_fraction, pfr_fast_recycle_fraction,
                        sol_fast_recycle_fraction, target_fast_recycle_energy_factor,
                        sol_fast_recycle_energy_factor, pfr_fast_recycle_energy_factor});
    // FIXME: These are global settings, but are being overwritten by each particular
    // recycling channel

    // Boolean flags for enabling recycling in different regions
    target_recycle = from_options["target_recycle"]
                         .doc("Recycling in the targets?")
                         .withDefault<bool>(false);

    sol_recycle = from_options["sol_recycle"]
                      .doc("Recycling in the SOL edge?")
                      .withDefault<bool>(false);

    pfr_recycle = from_options["pfr_recycle"]
                      .doc("Recycling in the PFR edge?")
                      .withDefault<bool>(false);

    neutral_pump = from_options["neutral_pump"]
                       .doc("Neutral pump enabled? Note, need location in grid file")
                       .withDefault<bool>(false);
  }

  if (target_recycle) {
    setPermissions(readIfSet("species:{from}:energy_flow_ylow"));
  }
  if (sol_recycle or pfr_recycle) {
    setPermissions(readIfSet("species:{from}:energy_flow_xlow"));
    setPermissions(readIfSet("species:{from}:particle_flow_xlow"));
  }
  substitutePermissions("to",
                        std::vector<std::string>(to_species.begin(), to_species.end()));
  substitutePermissions(
      "from", std::vector<std::string>(from_species.begin(), from_species.end()));
  substitutePermissions("to_inputs", {"AA", "density", "pressure", "temperature"});
  substitutePermissions("from_inputs", {"density", "velocity", "temperature"});
  substitutePermissions("outputs", {"density_source", "energy_source"});
}

void Recycling::transform_impl(GuardedOptions& state) {

  // Get metric tensor components
  Coordinates* coord = mesh->getCoordinates();
  const Field2D& J = coord->J;
  const Field2D& dy = coord->dy;
  const Field2D& dx = coord->dx;
  const Field2D& dz = coord->dz;
  const Field2D& g_22 = coord->g_22;

  for (auto& channel : channels) {
    const GuardedOptions species_from = state["species"][channel.from];

    const Field3D N = get<Field3D>(species_from["density"]);
    const Field3D V = get<Field3D>(species_from["velocity"]);    // Parallel flow velocity
    const Field3D T = get<Field3D>(species_from["temperature"]); // Ion temperature

    GuardedOptions species_to = state["species"][channel.to];
    const Field3D Nn = get<Field3D>(species_to["density"]);
    const Field3D Pn = get<Field3D>(species_to["pressure"]);
    const Field3D Tn = get<Field3D>(species_to["temperature"]);
    const BoutReal AAn = get<BoutReal>(species_to["AA"]);

    const Field3D Nnlim = floor(Nn, density_floor);
    const Field3D Pnlim = floor(Pn, pressure_floor);
    const Field3D Tnlim = Pnlim / Nnlim;

    // Recycling particle and energy sources will be added to these global sources
    // which are then passed to the density and pressure equations
    density_source = species_to.isSet("density_source")
                         ? getNonFinal<Field3D>(species_to["density_source"])
                         : 0.0;
    energy_source = species_to.isSet("energy_source")
                        ? getNonFinal<Field3D>(species_to["energy_source"])
                        : 0.0;
    debug = zeroFrom(Nn);

    // Recycling at the divertor target plates
    if (target_recycle) {

      if (!has_sheath_boundary_simple and (channel.target_fast_recycle_fraction > 0.0)) {
        throw BoutException(
            "Currently only sheath_boundary_simple is supported for fast recycling.");
      }

      if (has_sheath_boundary_simple) {
        // Fast recycling needs to know how much energy the "from" species is losing to
        // the boundary
        if (species_from.isSet("energy_flow_ylow")) {
          energy_flow_ylow = get<Field3D>(species_from["energy_flow_ylow"]);
        } else {
          energy_flow_ylow = 0.0;
        }
      }

      channel.target_recycle_density_source = 0.0;
      channel.target_recycle_energy_source = 0.0;

      // Lower Y boundary

      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {

          auto i = indexAt(Nn, r.ind, mesh->ystart, jz); // Final domain cell
          auto ig = i.ym();                              // First guard cell

          // Reproduce sheath advected particle flux from evolve_density
          BoutReal flux = -0.5 * (N[i] + N[ig]) * 0.5 * (V[i] + V[ig]);
          flux = std::max(flux, 0.0);

          BoutReal daparsheath = (J[i] + J[ig]) / (sqrt(g_22[i]) + sqrt(g_22[ig])) * 0.5
                                 * (dx[i] + dx[ig]) * 0.5 * (dz[i] + dz[ig]);
          BoutReal volume = J[i] * dx[i] * dy[i] * dz[i];

          // Flow of recycled neutrals into domain [s-1]
          BoutReal flow = channel.target_multiplier * flux * daparsheath;

          // Calculate sources in the final cell [m^-3 s^-1]
          channel.target_recycle_density_source[i] += flow / volume; // For diagnostic
          density_source[i] += flow / volume;                        // For use in solver

          // Reproduce sheath advected energy flux from evolve_pressure
          BoutReal nisheath = (N[i] + N[ig]) * 0.5;
          BoutReal tisheath = (T[i] + T[ig]) * 0.5;
          BoutReal visheath = (V[i] + V[ig]) * 0.5;
          BoutReal sheath_ion_heat_flow =
              abs(gamma_i * nisheath * tisheath * visheath * daparsheath / volume);

          // Blend fast (ion energy) and thermal (constant energy) recycling
          // Calculate returning neutral heat flow in [W]
          BoutReal recycle_energy_flow =
              // Fast recycling part
              sheath_ion_heat_flow * channel.target_multiplier
                  * channel.target_fast_recycle_energy_factor
                  * channel.target_fast_recycle_fraction
              // Thermal recycling part
              + flow * (1 - channel.target_fast_recycle_fraction) * channel.target_energy;

          // Divide heat flow in [W] by cell volume to get source in [m^-3 s^-1]
          channel.target_recycle_energy_source[i] += recycle_energy_flow / volume;
          energy_source[i] += recycle_energy_flow / volume;
        }
      }

      // Upper Y boundary

      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {

        for (int jz = 0; jz < mesh->LocalNz; jz++) {

          auto i = indexAt(Nn, r.ind, mesh->yend, jz); // Final domain cell
          auto ig = i.yp();                            // First guard cell

          // Reproduce sheath advected particle flux from evolve_density
          BoutReal flux = 0.5 * (N[i] + N[ig]) * 0.5 * (V[i] + V[ig]);
          flux = std::max(flux, 0.0);

          BoutReal daparsheath = (J[i] + J[ig]) / (sqrt(g_22[i]) + sqrt(g_22[ig])) * 0.5
                                 * (dx[i] + dx[ig]) * 0.5 * (dz[i] + dz[ig]);
          BoutReal volume = J[i] * dx[i] * dy[i] * dz[i];

          // Flow of recycled neutrals into domain [s-1]
          BoutReal flow = channel.target_multiplier * flux * daparsheath;

          // Calculate sources in the final cell [m^-3 s^-1]
          channel.target_recycle_density_source[i] += flow / volume; // For diagnostic
          density_source[i] += flow / volume;                        // For use in solver

          // Reproduce sheath advected energy flux from evolve_pressure
          BoutReal nisheath = (N[i] + N[ig]) * 0.5;
          BoutReal tisheath = (T[i] + T[ig]) * 0.5;
          BoutReal visheath = (V[i] + V[ig]) * 0.5;
          BoutReal sheath_ion_heat_flow =
              abs(gamma_i * nisheath * tisheath * visheath * daparsheath / volume);

          // Blend fast (ion energy) and thermal (constant energy) recycling
          // Calculate returning neutral heat flow in [W]
          BoutReal recycle_energy_flow =
              // Fast recycling part
              sheath_ion_heat_flow * channel.target_multiplier
                  * channel.target_fast_recycle_energy_factor
                  * channel.target_fast_recycle_fraction
              // Thermal recycling part
              + flow * (1 - channel.target_fast_recycle_fraction) * channel.target_energy;

          // Divide heat flow in [W] by cell volume to get source in [m^-3 s^-1]
          channel.target_recycle_energy_source[i] += recycle_energy_flow / volume;
          energy_source[i] += recycle_energy_flow / volume;
        }
      }
    }

    // Get edge particle and heat for the species being recycled
    if (sol_recycle or pfr_recycle) {

      // Initialise counters of pump recycling fluxes
      channel.wall_recycle_density_source = 0;
      channel.wall_recycle_energy_source = 0;
      channel.pump_density_source = 0;
      channel.pump_energy_source = 0;

      if (species_from.isSet("energy_flow_xlow")) {
        energy_flow_xlow = get<Field3D>(species_from["energy_flow_xlow"]);
      } else if ((channel.sol_fast_recycle_fraction > 0)
                 or (channel.pfr_fast_recycle_fraction > 0)) {
        throw BoutException("SOL/PFR fast recycle enabled but no cell edge heat flow "
                            "available, check your wall BC choice");
      };

      if (species_from.isSet("particle_flow_xlow")) {
        particle_flow_xlow = get<Field3D>(species_from["particle_flow_xlow"]);
      } else if ((channel.sol_fast_recycle_fraction > 0)
                 or (channel.pfr_fast_recycle_fraction > 0)) {
        throw BoutException("SOL/PFR fast recycle enabled but no cell edge particle flow "
                            "available, check your wall BC choice");
      };
    }

    // Recycling at the SOL edge (2D/3D only)
    if (sol_recycle) {

      // Flow out of domain is positive in the positive coordinate direction
      Field3D radial_particle_outflow = particle_flow_xlow;
      Field3D radial_energy_outflow = energy_flow_xlow;

      if (mesh->lastX()) { // Only do this for the processor which has the edge region
        for (int iy = 0; iy < mesh->LocalNy; iy++) {
          for (int iz = 0; iz < mesh->LocalNz; iz++) {

            // Volume of cell adjacent to wall which will receive source
            BoutReal volume = J(mesh->xend, iy) * dx(mesh->xend, iy) * dy(mesh->xend, iy)
                              * dz(mesh->xend, iy);

            // If cell is a pump, overwrite multiplier with pump multiplier
            BoutReal multiplier = channel.sol_multiplier;
            if ((is_pump(mesh->xend, iy) == 1.0) and (neutral_pump)) {
              multiplier = channel.pump_multiplier;
            }

            // Flow of recycled species back from the edge due to recycling
            // SOL edge = LHS flow of inner guard cells on the high X side
            // (mesh->xend+1) Recycling source is 0 for each cell where the flow goes
            // into instead of out of the domain
            BoutReal recycle_particle_flow = 0;
            if (radial_particle_outflow(mesh->xend + 1, iy, iz) > 0) {
              recycle_particle_flow =
                  multiplier * radial_particle_outflow(mesh->xend + 1, iy, iz);
            }

            BoutReal ion_energy_flow = radial_energy_outflow(
                mesh->xend + 1, iy, iz); // Ion heat flow to wall in [W]. This is on
                                         // xlow edge so take guard cell

            // Blend fast (ion energy) and thermal (constant energy) recycling
            // Calculate returning neutral heat flow in [W]
            BoutReal recycle_energy_flow =
                ion_energy_flow * channel.sol_multiplier
                    * channel.sol_fast_recycle_energy_factor
                    * channel.sol_fast_recycle_fraction // Fast recycling part
                + recycle_particle_flow * (1 - channel.sol_fast_recycle_fraction)
                      * channel.sol_energy; // Thermal recycling part

            // Calculate neutral pump neutral sinks due to neutral impingement
            BoutReal pump_neutral_particle_sink = 0.0;
            BoutReal pump_neutral_energy_sink = 0.0;

            // These are additional to particle sinks due to recycling
            if ((is_pump(mesh->xend, iy) == 1.0) and (neutral_pump)) {

              auto i = indexAt(Nn, mesh->xend, iy, iz); // Final domain cell
              auto is = i.xm();                         // Second to final domain cell
              auto ig = i.xp();                         // First guard cell

              // Free boundary condition on Nn, Pn, Tn
              // These are NOT communicated back into state and will exist only in this
              // component This will prevent neutrals leaking through cross-field
              // transport from neutral_mixed or other components While enabling us to
              // still calculate radial wall fluxes separately here
              BoutReal nnguard = SQ(Nn[i]) / Nnlim[is];
              BoutReal tnguard = SQ(Tn[i]) / Tnlim[is];

              // Calculate wall conditions
              BoutReal nnsheath = 0.5 * (Nn[i] + nnguard);
              BoutReal tnsheath = 0.5 * (Tn[i] + tnguard);
              BoutReal v_th =
                  0.25
                  * sqrt(8 * tnsheath / (PI * AAn)); // Stangeby p.69 eqns. 2.21, 2.24

              // Convert dy to poloidal length: dl = dy * sqrt(g22) = dy * h_theta
              // Convert dz to toroidal length:  = dz*sqrt(g_33) = dz * R = 2piR
              // Calculate radial wall area in [m^2]
              // Calculate final cell volume [m^3]
              BoutReal dpolsheath =
                  0.5 * (coord->dy[i] + coord->dy[ig]) * 1
                  / (0.5 * (sqrt(coord->g22[i]) + sqrt(coord->g22[ig])));
              BoutReal dtorsheath = 0.5 * (coord->dz[i] + coord->dz[ig]) * 0.5
                                    * (sqrt(coord->g_33[i]) + sqrt(coord->g_33[ig]));
              BoutReal dasheath = dpolsheath * dtorsheath; // [m^2]
              BoutReal dv = coord->J[i] * coord->dx[i] * coord->dy[i] * coord->dz[i];

              // Calculate particle and energy fluxes of neutrals hitting the pump
              // Assume thermal velocity greater than perpendicular velocity and use it
              // for flux calc
              BoutReal pump_neutral_particle_flow = v_th * nnsheath * dasheath; // [s^-1]
              pump_neutral_particle_sink =
                  pump_neutral_particle_flow / dv
                  * (1 - multiplier); // Particle sink [s^-1 m^-3]

              // Use gamma=2.0 as per Stangeby p.69, total energy of static Maxwellian
              BoutReal pump_neutral_energy_flow =
                  2.0 * tnsheath * v_th * nnsheath * dasheath; // [W]
              pump_neutral_energy_sink =
                  pump_neutral_energy_flow / dv * (1 - multiplier); // heatsink [W m^-3]

              // Divide flows by volume to get sources
              // Save to pump diagnostic
              channel.pump_density_source(mesh->xend, iy, iz) +=
                  recycle_particle_flow / volume - pump_neutral_particle_sink;
              channel.pump_energy_source(mesh->xend, iy, iz) +=
                  recycle_energy_flow / volume - pump_neutral_energy_sink;

            } else {
              // Save to wall diagnostic (pump sinks are 0 if not on pump)
              channel.wall_recycle_density_source(mesh->xend, iy, iz) +=
                  recycle_particle_flow / volume;
              channel.wall_recycle_energy_source(mesh->xend, iy, iz) +=
                  recycle_energy_flow / volume;
            }

            // Save to solver
            density_source(mesh->xend, iy, iz) +=
                recycle_particle_flow / volume - pump_neutral_particle_sink;
            energy_source(mesh->xend, iy, iz) +=
                recycle_energy_flow / volume - pump_neutral_energy_sink;
          }
        }
      }
    }

    // Recycling at the PFR edge (2D/3D only)
    if (pfr_recycle) {

      // PFR is flipped compared to edge: x=0 is at the PFR edge. Therefore outflow is
      // in the negative coordinate direction.
      Field3D radial_particle_outflow = particle_flow_xlow * -1;
      Field3D radial_energy_outflow = energy_flow_xlow * -1;

      if (mesh->firstX()) { // Only do this for the processor which has the core region
        if (!mesh->periodicY(mesh->xstart)) { // Only do this for the processor with a
                                              // periodic Y, i.e. the PFR
          for (int iy = 0; iy < mesh->LocalNy; iy++) {
            for (int iz = 0; iz < mesh->LocalNz; iz++) {

              // Volume of cell adjacent to wall which will receive source
              BoutReal volume = J(mesh->xstart, iy) * dx(mesh->xstart, iy)
                                * dy(mesh->xstart, iy) * dz(mesh->xstart, iy);

              // If cell is a pump, overwrite multiplier with pump multiplier
              BoutReal multiplier = channel.pfr_multiplier;
              if ((is_pump(mesh->xstart, iy) == 1.0) and (neutral_pump)) {
                multiplier = channel.pump_multiplier;
              }

              // Flow of recycled species back from the edge due to recycling
              // PFR edge = LHS flow of the first domain cell on the low X side
              // (mesh->xstart) Recycling source is 0 for each cell where the flow goes
              // into instead of out of the domain
              BoutReal recycle_particle_flow = 0;
              if (radial_particle_outflow(mesh->xstart, iy, iz) > 0) {
                recycle_particle_flow =
                    multiplier * radial_particle_outflow(mesh->xstart, iy, iz);
              }

              BoutReal ion_energy_flow = radial_energy_outflow(
                  mesh->xstart, iy, iz); // Ion heat flow to wall in [W]. This is on
                                         // xlow edge so take first domain cell

              // Blend fast (ion energy) and thermal (constant energy) recycling
              // Calculate returning neutral heat flow in [W]
              BoutReal recycle_energy_flow =
                  ion_energy_flow * channel.pfr_multiplier
                      * channel.pfr_fast_recycle_energy_factor
                      * channel.pfr_fast_recycle_fraction // Fast recycling part
                  + recycle_particle_flow * (1 - channel.pfr_fast_recycle_fraction)
                        * channel.pfr_energy; // Thermal recycling part

              // Calculate neutral pump neutral sinks due to neutral impingement
              // These are additional to particle sinks due to recycling
              BoutReal pump_neutral_particle_sink = 0.0;
              BoutReal pump_neutral_energy_sink = 0.0;

              // Add to appropriate diagnostic field depending if pump or not
              if ((is_pump(mesh->xstart, iy) == 1.0) and (neutral_pump)) {

                auto i = indexAt(Nn, mesh->xstart, iy, iz); // Final domain cell
                auto is = i.xp();                           // Second to final domain cell
                auto ig = i.xm();                           // First guard cell

                // Free boundary condition on Nn, Pn, Tn
                // These are NOT communicated back into state and will exist only in
                // this component This will prevent neutrals leaking through cross-field
                // transport from neutral_mixed or other components While enabling us to
                // still calculate radial wall fluxes separately here
                BoutReal nnguard = SQ(Nn[i]) / Nnlim[is];
                BoutReal tnguard = SQ(Tn[i]) / Tnlim[is];

                // Calculate wall conditions
                BoutReal nnsheath = 0.5 * (Nn[i] + nnguard);
                BoutReal tnsheath = 0.5 * (Tn[i] + tnguard);
                BoutReal v_th =
                    0.25
                    * sqrt(8 * tnsheath / (PI * AAn)); // Stangeby p.69 eqns. 2.21, 2.24

                // Convert dy to poloidal length: dl = dy * sqrt(g22) = dy * h_theta
                // Convert dz to toroidal length:  = dz*sqrt(g_33) = dz * R = 2piR
                // Calculate radial wall area in [m^2]
                // Calculate final cell volume [m^3]
                BoutReal dpolsheath =
                    0.5 * (coord->dy[i] + coord->dy[ig]) * 1
                    / (0.5 * (sqrt(coord->g22[i]) + sqrt(coord->g22[ig])));
                BoutReal dtorsheath = 0.5 * (coord->dz[i] + coord->dz[ig]) * 0.5
                                      * (sqrt(coord->g_33[i]) + sqrt(coord->g_33[ig]));
                BoutReal dasheath = dpolsheath * dtorsheath; // [m^2]
                BoutReal dv = coord->J[i] * coord->dx[i] * coord->dy[i] * coord->dz[i];

                // Calculate particle and energy fluxes of neutrals hitting the pump
                // Assume thermal velocity greater than perpendicular velocity and use
                // it for flux calc
                BoutReal pump_neutral_particle_flow =
                    v_th * nnsheath * dasheath; // [s^-1]
                pump_neutral_particle_sink =
                    pump_neutral_particle_flow / dv
                    * (1 - multiplier); // Particle sink [s^-1 m^-3]

                // Use gamma=2.0 as per Stangeby p.69, total energy of static Maxwellian
                BoutReal pump_neutral_energy_flow =
                    2.0 * tnsheath * v_th * nnsheath * dasheath; // [W]
                pump_neutral_energy_sink =
                    pump_neutral_energy_flow / dv * (1 - multiplier); // heatsink [W m^-3]

                // Divide flows by volume to get sources
                // Save to pump diagnostic
                channel.pump_density_source(mesh->xstart, iy, iz) +=
                    recycle_particle_flow / volume - pump_neutral_particle_sink;
                channel.pump_energy_source(mesh->xstart, iy, iz) +=
                    recycle_energy_flow / volume - pump_neutral_energy_sink;

              } else {
                // Save to wall diagnostic (pump sinks are 0 if not on pump)
                channel.wall_recycle_density_source(mesh->xstart, iy, iz) +=
                    recycle_particle_flow / volume - pump_neutral_particle_sink;
                channel.wall_recycle_energy_source(mesh->xstart, iy, iz) +=
                    recycle_energy_flow / volume - pump_neutral_energy_sink;
              }

              // Save to solver
              density_source(mesh->xstart, iy, iz) +=
                  recycle_particle_flow / volume - pump_neutral_particle_sink;
              energy_source(mesh->xstart, iy, iz) +=
                  recycle_energy_flow / volume - pump_neutral_energy_sink;
            }
          }
        }
      }
    }

    // Put the updated sources back into the state
    set<Field3D>(species_to["density_source"], density_source);
    set<Field3D>(species_to["energy_source"], energy_source);
  }
}

void Recycling::outputVars(Options& state) {

  if (neutral_pump) {
    // Save the pump mask as a time-independent field
    set_with_attrs(state[{std::string("is_pump")}], is_pump,
                   {{"standard_name", "neutral pump location"},
                    {"long_name", std::string("Neutral pump location")},
                    {"source", "recycling"}});
  }

  if (diagnose) {
    // Normalisations
    auto Nnorm = get<BoutReal>(state["Nnorm"]);
    auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
    auto Tnorm = get<BoutReal>(state["Tnorm"]);
    BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

    for (const auto& channel : channels) {
      // Save particle and energy source for the species created during recycling
      set_with_attrs(state["debug"], debug,
                     {{"time_dimension", "t"},
                      {"units", "-"},
                      {"conversion", 1},
                      {"standard_name", "debug variable"},
                      {"long_name", "debug variable"},
                      {"source", "recycling"}});
      // Target recycling
      if (target_recycle) {
        set_with_attrs(
            state[{std::string("S") + channel.to + std::string("_target_recycle")}],
            channel.target_recycle_density_source,
            {{"time_dimension", "t"},
             {"units", "m^-3 s^-1"},
             {"conversion", Nnorm * Omega_ci},
             {"standard_name", "particle source"},
             {"long_name",
              std::string("Target recycling particle source of ") + channel.to},
             {"source", "recycling"}});

        set_with_attrs(
            state[{std::string("E") + channel.to + std::string("_target_recycle")}],
            channel.target_recycle_energy_source,
            {{"time_dimension", "t"},
             {"units", "W m^-3"},
             {"conversion", Pnorm * Omega_ci},
             {"standard_name", "energy source"},
             {"long_name",
              std::string("Target recycling energy source of ") + channel.to},
             {"source", "recycling"}});
      }

      // Wall recycling
      if ((sol_recycle) or (pfr_recycle)) {
        set_with_attrs(
            state[{std::string("S") + channel.to + std::string("_wall_recycle")}],
            channel.wall_recycle_density_source,
            {{"time_dimension", "t"},
             {"units", "m^-3 s^-1"},
             {"conversion", Nnorm * Omega_ci},
             {"standard_name", "particle source"},
             {"long_name",
              std::string("Wall recycling particle source of ") + channel.to},
             {"source", "recycling"}});

        set_with_attrs(
            state[{std::string("E") + channel.to + std::string("_wall_recycle")}],
            channel.wall_recycle_energy_source,
            {{"time_dimension", "t"},
             {"units", "W m^-3"},
             {"conversion", Pnorm * Omega_ci},
             {"standard_name", "energy source"},
             {"long_name", std::string("Wall recycling energy source of ") + channel.to},
             {"source", "recycling"}});
      }

      // Neutral pump
      if (neutral_pump) {

        if (channel.pump_density_source.isAllocated()) {
          set_with_attrs(state[{std::string("S") + channel.to + std::string("_pump")}],
                         channel.pump_density_source,
                         {{"time_dimension", "t"},
                          {"units", "m^-3 s^-1"},
                          {"conversion", Nnorm * Omega_ci},
                          {"standard_name", "particle source"},
                          {"long_name", std::string("Pump recycling particle source of ")
                                            + channel.to},
                          {"source", "recycling"}});

          set_with_attrs(state[{std::string("E") + channel.to + std::string("_pump")}],
                         channel.pump_energy_source,
                         {{"time_dimension", "t"},
                          {"units", "W m^-3"},
                          {"conversion", Pnorm * Omega_ci},
                          {"standard_name", "energy source"},
                          {"long_name",
                           std::string("Pump recycling energy source of ") + channel.to},
                          {"source", "recycling"}});
        } else {
          throw BoutException(
              "Error: neutral pump sources unallocated, likely because recycling "
              "on the relevant surface is disabled. Check that all boundaries with a "
              "neutral pump "
              "have recycling enabled!");
        }
      }
    }
  }
}
