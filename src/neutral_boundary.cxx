#include "hermes_utils.hxx"
#include "bout/mesh.hxx"
#include <bout/constants.hxx>
using bout::globals::mesh;

#include "../include/neutral_boundary.hxx"

NeutralBoundary::NeutralBoundary(std::string name, Options& alloptions,
                                 [[maybe_unused]] Solver* solver)
    : Component({writeBoundary("species:{name}:{boundary_outputs}"),
                 writeBoundaryIfSet("species:{name}:{conditional_outputs}"),
                 readWrite("species:{name}:energy_source"),
                 readOnly("species:{name}:AA"),
                 readOnly("species:{from}:{from_inputs}"),
                 readOnly("species:{to}:{to_inputs}"),
                 readWrite("species:{to}:{outputs}"),
                 readWrite("species:{from}:{outputs}")}),
      name(name) {

  auto& options = alloptions[name];
  const Options& units = alloptions["units"];
  Tnorm = units["eV"];

  diagnose =
      options["diagnose"].doc("Save additional diagnostics?").withDefault<bool>(false);
  lower_y = options["neutral_boundary_lower_y"]
                .doc("Boundary on lower y?")
                .withDefault<bool>(true);
  upper_y = options["neutral_boundary_upper_y"]
                .doc("Boundary on upper y?")
                .withDefault<bool>(true);
  sol = options["neutral_boundary_sol"].doc("Boundary on SOL?").withDefault<bool>(false);
  pfr = options["neutral_boundary_pfr"].doc("Boundary on PFR?").withDefault<bool>(false);

  target_energy_refl_factor = options["target_energy_refl_factor"]
                                  .doc("Fraction of energy retained by neutral particles "
                                       "after wall reflection at target")
                                  .withDefault<BoutReal>(0.75);

  sol_energy_refl_factor = options["sol_energy_refl_factor"]
                               .doc("Fraction of energy retained by neutral particles "
                                    "after wall reflection at SOL")
                               .withDefault<BoutReal>(0.75);

  pfr_energy_refl_factor = options["pfr_energy_refl_factor"]
                               .doc("Fraction of energy retained by neutral particles "
                                    "after wall reflection at PFR")
                               .withDefault<BoutReal>(0.75);

  target_fast_refl_fraction =
      options["target_fast_refl_fraction"]
          .doc("Fraction of neutrals that are undergoing fast reflection at the target")
          .withDefault<BoutReal>(0.8);

  sol_fast_refl_fraction =
      options["sol_fast_refl_fraction"]
          .doc("Fraction of neutrals that are undergoing fast reflection at the sol")
          .withDefault<BoutReal>(0.8);

  pfr_fast_refl_fraction =
      options["pfr_fast_refl_fraction"]
          .doc("Fraction of neutrals that are undergoing fast reflection at the pfr")
          .withDefault<BoutReal>(0.8);

  /// For neutral ionising at core
  std::set<std::string> from_species, to_species;

  density_floor =
      options["density_floor"].doc("Minimum density floor").withDefault(1e-7);
  pressure_floor = density_floor * (1. / get<BoutReal>(alloptions["units"]["eV"]));

  core_ionising = options["core_ionising"]
                    .doc("Neutrals ionised in the core?")
                    .withDefault<bool>(false);

  if (core_ionising) {
    std::string from = name;
    std::string to = options["ionise_as"]
                         .doc("Name of the species to ionise into")
                         .as<std::string>();

    from_species.insert(from);
    to_species.insert(to);

    BoutReal core_ionise_multiplier =
        options["core_ionise_multiplier"]
            .doc("Multiply the core ionised flux by this factor. Should be >=0 and <= 1")
            .withDefault<BoutReal>(1.0);
    // BoutReal core_ionise_energy =
    //     options["core_ionise_energy"]
    //         .doc("Fixed energy of the ionised particles at core [eV]")
    //         .withDefault<BoutReal>(3.0)
    //     / Tnorm; // Normalise from eV
    if ((core_ionise_multiplier < 0.0) or (core_ionise_multiplier > 1.0)){
      throw BoutException("Core ionise multipliers must be betweeen 0 and 1");
    }
    channels.push_back({from, to, core_ionise_multiplier, 0.0});
  }
  
  substitutePermissions("name", {name});
  substitutePermissions("boundary_outputs", {"density", "temperature", "pressure"});
  substitutePermissions("conditional_outputs", {"velocity", "momentum"});
  // For core ionising
  if (core_ionising) {
    setPermissions(readIfSet("species:{from}:energy_flow_xlow"));
    setPermissions(readIfSet("species:{from}:particle_flow_xlow"));
  }
  substitutePermissions("to",
                        std::vector<std::string>(to_species.begin(), to_species.end()));
  substitutePermissions(
      "from", std::vector<std::string>(from_species.begin(), from_species.end()));
  substitutePermissions("to_inputs", {"AA", "density", "pressure", "temperature"});
  substitutePermissions("from_inputs", {"density", "velocity", "temperature"});
  substitutePermissions("outputs", {"density_source", "momentum_source", "energy_source"});
}


void NeutralBoundary::transform_impl(GuardedOptions& state) {
  auto species = state["species"][name];
  const BoutReal AA = get<BoutReal>(species["AA"]);

  Field3D Nn = toFieldAligned(GET_NOBOUNDARY(Field3D, species["density"]));
  Field3D Pn = toFieldAligned(GET_NOBOUNDARY(Field3D, species["pressure"]));
  Field3D Tn = toFieldAligned(GET_NOBOUNDARY(Field3D, species["temperature"]));

  Field3D Vn = IS_SET_NOBOUNDARY(species["velocity"])
                   ? toFieldAligned(getNoBoundary<Field3D>(species["velocity"]))
                   : zeroFrom(Nn);

  Field3D NVn = IS_SET_NOBOUNDARY(species["momentum"])
                    ? toFieldAligned(getNoBoundary<Field3D>(species["momentum"]))
                    : zeroFrom(Nn);

  // Get the energy source, or create if not set
  Field3D energy_source =
      species.isSet("energy_source")
          ? toFieldAligned(getNonFinal<Field3D>(species["energy_source"]))
          : zeroFrom(Nn);

  Coordinates* coord = mesh->getCoordinates();
  const Field2D& J = coord->J;
  const Field2D& dy = coord->dy;
  const Field2D& dx = coord->dx;
  const Field2D& dz = coord->dz;
  target_energy_source = 0;
  wall_energy_source = 0;

  // Targets
  if (lower_y) {
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(Nn, r.ind, mesh->ystart, jz);
        auto im = i.ym();

        // Free boundary condition on Nn, Pn, Tn
        // This is problematic when Nn, Pn or Tn are zero
        // Nn[im] = SQ(Nn[i]) / Nn[ip];
        // Pn[im] = SQ(Pn[i]) / Pn[ip];
        // Tn[im] = SQ(Tn[i]) / Tn[ip];

        // Neumann boundary condition: do not extrapolate, but
        // assume the target value is same as the final cell centre.
        // Shouldn't affect results much and more resilient to positivity issues
        Nn[im] = Nn[i];
        Pn[im] = Pn[i];
        Tn[im] = Tn[i];

        // No-flow boundary condition
        Vn[im] = -Vn[i];
        NVn[im] = -NVn[i];

        // Calculate midpoint values at wall
        const BoutReal nnsheath = 0.5 * (Nn[im] + Nn[i]);
        const BoutReal tnsheath = 0.5 * (Tn[im] + Tn[i]);

        // Thermal speed
        const BoutReal v_th =
            0.25 * sqrt(8 * tnsheath / (PI * AA)); // Stangeby p.69 eqns. 2.21, 2.24

        // Approach adapted from D. Power thesis 2023
        BoutReal T_FC = 3 / Tnorm; // Franck-Condon temp (hardcoded for now)

        // Outgoing neutral heat flux [W/m^2]
        // This is rearranged from Power for clarity - note definition of v_th.
        // Uses standard Stangeby 1D static Maxwellian particle/heat fluxes for fast terms
        // and simply Q = T * particle flux for the monoenergetic thermal reflected
        // population.
        BoutReal q = 2 * nnsheath * tnsheath * v_th // Incident energy
                     - (target_energy_refl_factor * target_fast_refl_fraction) * 2
                           * nnsheath * tnsheath * v_th // Fast reflected energy
                     - (1 - target_fast_refl_fraction) * T_FC * nnsheath
                           * v_th; // Thermal reflected energy

        // Cross-sectional area in XZ plane:
        BoutReal da = (coord->J[i] + coord->J[im])
                      / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[im])) * 0.5
                      * (coord->dx[i] + coord->dx[im]) * 0.5
                      * (coord->dz[i] + coord->dz[im]); // [m^2]

        // Multiply by area to get energy flow (power)
        BoutReal flow = q * da; // [W]

        // Divide by cell volume to get source [W/m^3]
        BoutReal cooling_source =
            flow / (coord->dx[i] * coord->dy[i] * coord->dz[i] * coord->J[i]);

        // Subtract from cell next to boundary
        energy_source[i] -= cooling_source;
        target_energy_source[i] -= cooling_source;
      }
    }
  }

  if (upper_y) {
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        auto i = indexAt(Nn, r.ind, mesh->yend, jz);
        auto ip = i.yp();

        // Free boundary condition on Nn, Pn, Tn
        // This is problematic when Nn, Pn or Tn are zero
        // Nn[ip] = SQ(Nn[i]) / Nn[im];
        // Pn[ip] = SQ(Pn[i]) / Pn[im];
        // Tn[ip] = SQ(Tn[i]) / Tn[im];

        // Neumann boundary condition: do not extrapolate, but
        // assume the target value is same as the final cell centre.
        // Shouldn't affect results much and more resilient to positivity issues
        Nn[ip] = Nn[i];
        Pn[ip] = Pn[i];
        Tn[ip] = Tn[i];

        // No-flow boundary condition
        Vn[ip] = -Vn[i];
        NVn[ip] = -NVn[i];

        // Calculate midpoint values at wall
        const BoutReal nnsheath = 0.5 * (Nn[ip] + Nn[i]);
        const BoutReal tnsheath = 0.5 * (Tn[ip] + Tn[i]);

        // Thermal speed
        const BoutReal v_th =
            0.25 * sqrt(8 * tnsheath / (PI * AA)); // Stangeby p.69 eqns. 2.21, 2.24

        // Approach adapted from D. Power thesis 2023
        BoutReal T_FC = 3 / Tnorm; // Franck-Condon temp (hardcoded for now)

        // Outgoing neutral heat flux [W/m^2]
        // This is rearranged from Power for clarity - note definition of v_th.
        BoutReal q = 2 * nnsheath * tnsheath * v_th // Incident energy
                     - (target_energy_refl_factor * target_fast_refl_fraction) * 2
                           * nnsheath * tnsheath * v_th // Fast reflected energy
                     - (1 - target_fast_refl_fraction) * T_FC * nnsheath
                           * v_th; // Thermal reflected energy

        // Cross-sectional area in XZ plane:
        BoutReal da = (coord->J[i] + coord->J[ip])
                      / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[ip])) * 0.5
                      * (coord->dx[i] + coord->dx[ip]) * 0.5
                      * (coord->dz[i] + coord->dz[ip]); // [m^2]

        // Multiply by area to get energy flow (power)
        BoutReal flow = q * da; // [W]

        // Divide by cell volume to get source [W/m^3]
        BoutReal cooling_source =
            flow / (coord->dx[i] * coord->dy[i] * coord->dz[i] * coord->J[i]);

        // Subtract from cell next to boundary
        energy_source[i] -= cooling_source;
        target_energy_source[i] -= cooling_source;
      }
    }
  }

  // SOL edge
  if (sol) {
    if (mesh->lastX()) { // Only do this for the processor which has the edge region
      for (int iy = 0; iy < mesh->LocalNy; iy++) {
        for (int iz = 0; iz < mesh->LocalNz; iz++) {

          auto i = indexAt(Nn, mesh->xend, iy, iz);      // Final domain cell
          auto ig = indexAt(Nn, mesh->xend + 1, iy, iz); // Guard cell

          // Calculate midpoint values at wall
          const BoutReal nnsheath = 0.5 * (Nn[ig] + Nn[i]);
          const BoutReal tnsheath = 0.5 * (Tn[ig] + Tn[i]);

          // Thermal speed of static Maxwellian in one direction
          const BoutReal v_th =
              0.25 * sqrt(8 * tnsheath / (PI * AA)); // Stangeby p.69 eqns. 2.21, 2.24

          // Approach adapted from D. Power thesis 2023
          BoutReal T_FC = 3 / Tnorm; // Franck-Condon temp (hardcoded for now)

          // Outgoing neutral heat flux [W/m^2]
          // This is rearranged from Power for clarity - note definition of v_th.
          BoutReal q = 2 * nnsheath * tnsheath * v_th // Incident energy
                       - (sol_energy_refl_factor * sol_fast_refl_fraction) * 2 * nnsheath
                             * tnsheath * v_th // Fast reflected energy
                       - (1 - sol_fast_refl_fraction) * T_FC * nnsheath
                             * v_th; // Thermal reflected energy

          // Multiply by radial cell area to get power
          // Expanded form of the calculation for clarity

          // Converts dy to poloidal length: dl = dy / sqrt(g22) = dy * h_theta
          BoutReal dpol = 0.5 * (coord->dy[i] + coord->dy[ig]) * 1
                          / (0.5 * (sqrt(coord->g22[i]) + sqrt(coord->g22[ig])));

          // Converts dz to toroidal length:  = dz*sqrt(g_33) = dz * R = 2piR
          BoutReal dtor = 0.5 * (coord->dz[i] + coord->dz[ig]) * 0.5
                          * (sqrt(coord->g_33[i]) + sqrt(coord->g_33[ig]));

          BoutReal da = dpol * dtor; // [m^2]

          // Multiply by area to get energy flow (power)
          BoutReal flow = q * da; // [W]

          // Divide by cell volume to get source [W/m^3]
          BoutReal cooling_source =
              flow
              / (coord->J[i] * coord->dx[i] * coord->dy[i] * coord->dz[i]); // [W m^-3]

          // Subtract from cell next to boundary
          energy_source[i] -= cooling_source;
          wall_energy_source[i] -= cooling_source;
        }
      }
    }
  }

  // PFR edge
  if (pfr) {
    if ((mesh->firstX())
        and (!mesh->periodicY(
            mesh->xstart))) { // do loop if inner edge and not periodic (i.e. PFR)
      for (int iy = 0; iy < mesh->LocalNy; iy++) {
        for (int iz = 0; iz < mesh->LocalNz; iz++) {

          auto i = indexAt(Nn, mesh->xstart, iy, iz);      // Final domain cell
          auto ig = indexAt(Nn, mesh->xstart - 1, iy, iz); // Guard cell

          // Calculate midpoint values at wall
          const BoutReal nnsheath = 0.5 * (Nn[ig] + Nn[i]);
          const BoutReal tnsheath = 0.5 * (Tn[ig] + Tn[i]);

          // Thermal speed of static Maxwellian in one direction
          const BoutReal v_th =
              0.25 * sqrt(8 * tnsheath / (PI * AA)); // Stangeby p.69 eqns. 2.21, 2.24

          // Approach adapted from D. Power thesis 2023
          BoutReal T_FC = 3 / Tnorm; // Franck-Condon temp (hardcoded for now)

          // Outgoing neutral heat flux [W/m^2]
          // This is rearranged from Power for clarity - note definition of v_th.
          BoutReal q = 2 * nnsheath * tnsheath * v_th // Incident energy
                       - (pfr_energy_refl_factor * pfr_fast_refl_fraction) * 2 * nnsheath
                             * tnsheath * v_th // Fast reflected energy
                       - (1 - pfr_fast_refl_fraction) * T_FC * nnsheath
                             * v_th; // Thermal reflected energy

          // Multiply by radial cell area to get power
          // Expanded form of the calculation for clarity

          // Converts dy to poloidal length: dl = dy / sqrt(g22) = dy * h_theta
          BoutReal dpol = 0.5 * (coord->dy[i] + coord->dy[ig]) * 1
                          / (0.5 * (sqrt(coord->g22[i]) + sqrt(coord->g22[ig])));

          // Converts dz to toroidal length:  = dz*sqrt(g_33) = dz * R = 2piR
          BoutReal dtor = 0.5 * (coord->dz[i] + coord->dz[ig]) * 0.5
                          * (sqrt(coord->g_33[i]) + sqrt(coord->g_33[ig]));

          BoutReal da = dpol * dtor; // [m^2]

          // Multiply by area to get energy flow (power)
          BoutReal flow = q * da; // [W]

          // Divide by cell volume to get source [W/m^3]
          BoutReal cooling_source =
              flow
              / (coord->J[i] * coord->dx[i] * coord->dy[i] * coord->dz[i]); // [W m^-3]

          // Subtract from cell next to boundary
          energy_source[i] -= cooling_source;
          wall_energy_source[i] -= cooling_source;
        }
      }
    }
  }

  // Set density, pressure and temperature, now with boundary conditions
  setBoundary(species["density"], fromFieldAligned(Nn));
  setBoundary(species["temperature"], fromFieldAligned(Tn));
  setBoundary(species["pressure"], fromFieldAligned(Pn));
  if (IS_SET_NOBOUNDARY(species["velocity"])) {
    setBoundary(species["velocity"], fromFieldAligned(Vn));
  }
  if (IS_SET_NOBOUNDARY(species["momentum"])) {
    setBoundary(species["momentum"], fromFieldAligned(NVn));
  }

  // Set energy source (negative in cell next to sheath)
  // Note: energy_source includes any sources previously set in other components
  set(species["energy_source"], fromFieldAligned(energy_source));



  /* ---- Core ionising ----*/
  if (core_ionising){ 
    for (auto& channel : channels) {

      GuardedOptions species_from = state["species"][channel.from];
      const Field3D Nn = get<Field3D>(species_from["density"]);
      const Field3D Pn = get<Field3D>(species_from["pressure"]);
      const Field3D Tn = get<Field3D>(species_from["temperature"]);
      const BoutReal AAn = get<BoutReal>(species_from["AA"]);
 
      const Field3D Nnlim = floor(Nn, density_floor);
      const Field3D Pnlim = floor(Pn, pressure_floor);
      const Field3D Tnlim = Pnlim / Nnlim;

      const GuardedOptions species_to = state["species"][channel.to];
      const Field3D N = get<Field3D>(species_to["density"]);
      const Field3D V = get<Field3D>(species_to["velocity"]);    // Parallel flow velocity
      const Field3D T = get<Field3D>(species_to["temperature"]); // Ion temperature

      // Recycling particle and energy sources will be added to these global sources
      // which are then passed to the density and pressure equations
      ion_density_source = species_to.isSet("density_source")
                          ? getNonFinal<Field3D>(species_to["density_source"])
                          : 0.0;
      ion_energy_source = species_to.isSet("energy_source")
                          ? getNonFinal<Field3D>(species_to["energy_source"])
                          : 0.0;
      ion_momentum_source = species_to.isSet("momentum_source")
                          ? getNonFinal<Field3D>(species_to["momentum_source"])
                          : 0.0;
      neutral_density_source = species_from.isSet("density_source")
                          ? getNonFinal<Field3D>(species_from["density_source"])
                          : 0.0;
      neutral_energy_source = species_from.isSet("energy_source")
                          ? getNonFinal<Field3D>(species_from["energy_source"])
                          : 0.0;
      neutral_momentum_source = species_from.isSet("momentum_source")
                          ? getNonFinal<Field3D>(species_from["momentum_source"])
                          : 0.0;
    
      channel.core_ion_density_source = 0.0;
      channel.core_ion_energy_source = 0.0;
      channel.core_ion_momentum_source = 0.0;
      channel.core_neutral_density_sink = 0.0;
      channel.core_neutral_energy_sink = 0.0;
      channel.core_neutral_momentum_sink = 0.0;
    
      if (species_from.isSet("energy_flow_xlow")) {
        energy_flow_xlow = get<Field3D>(species_from["energy_flow_xlow"]);
      } else {
        throw BoutException("core ionise enabled but no cell edge heat flow "
                            "available, check your core BC choice (?)");
      };
    
      if (species_from.isSet("particle_flow_xlow")) {
        particle_flow_xlow = get<Field3D>(species_from["particle_flow_xlow"]);
      } else {
        throw BoutException("core ionise enabled but no cell edge particle flow "
                            "available, check your core BC choice");
      };
    
      // Field3D radial_particle_outflow = particle_flow_xlow * -1.0;
      // Field3D radial_energy_outflow = energy_flow_xlow * -1.0;

      // The neutrals that reach the core boundary are removed from the domain 
      // and their flux is added to the boundary condition on the flux of ions
      //  from the core.
      // firstX from first domain cell 
      // the processor that contains the core boundary -> firstX
      // xstart -> the first domain cell in the boundary
      if (mesh->firstX() && mesh->periodicY(mesh->xstart)){
        for(int iy=0; iy < mesh->LocalNy; iy++){ 
          for(int iz = 0; iz < mesh->LocalNz; iz++){
            // Volume of cell adjacent to wall which will receive source
            BoutReal volume = J(mesh->xstart, iy) * dx(mesh->xstart, iy) * dy(mesh->xstart, iy)
                              * dz(mesh->xstart, iy);
            
            BoutReal multiplier = channel.core_multiplier;
            
            // Particle
            // Rather than using particle flow calculated in neutral transport,
            // use thermal speed to calculate particle flux into the core from 
            // core neutrals. This assumes that particles move thermally to 
            // ward the wall until they hit it.
            auto i = indexAt(Nn, mesh->xstart, iy, iz);      // Final domain cell
            auto ig = indexAt(Nn, mesh->xstart - 1, iy, iz); // Guard cell

            // Calculate midpoint values at wall
            const BoutReal nncore = 0.5 * (Nn[ig] + Nn[i]);
            const BoutReal tncore = 0.5 * (Tn[ig] + Tn[i]);

            // Thermal speed of static Maxwellian in one direction
            const BoutReal v_th =
                0.25 * sqrt(8 * tncore / (PI * AAn)); // Stangeby p.69 eqns. 2.21, 2.24
            
            // Convert dy to poloidal length: dl = dy * sqrt(g22) = dy * h_theta
            // Convert dz to toroidal length:  = dz*sqrt(g_33) = dz * R = 2piR
            // Calculate radial wall area in [m^2]
            // Calculate final cell volume [m^3]
            BoutReal dpolcore =
                0.5 * (coord->dy[i] + coord->dy[ig]) * 1
                / (0.5 * (sqrt(coord->g22[i]) + sqrt(coord->g22[ig])));
            BoutReal dtorcore = 0.5 * (coord->dz[i] + coord->dz[ig]) * 0.5
                                  * (sqrt(coord->g_33[i]) + sqrt(coord->g_33[ig]));
            BoutReal dacore = dpolcore * dtorcore; // [m^2]
            // BoutReal dv = coord->J[i] * coord->dx[i] * coord->dy[i] * coord->dz[i];
            BoutReal neutral_particle_flow_to_core = v_th * nncore * dacore;
            // The amount of neutral ionised in core
            BoutReal ionise_particle_flow = neutral_particle_flow_to_core * multiplier;
            
            // if (radial_particle_outflow(mesh->xstart-1, iy, iz) > 0.0) {
            //   ionise_particle_flow = 
            //     multiplier * radial_particle_outflow(mesh->xstart, iy, iz);
            // }
            // BoutReal core_neutral_particle_sink_flow = ionise_particle_flow;

            // Momentum 
            // mnv = kg m/s 1/m^3 
            BoutReal neutral_momentum_flow_to_core = v_th * v_th * nncore * dacore * AAn;
            BoutReal ionise_momentum_flow = neutral_momentum_flow_to_core * multiplier;

            // Energy 
            BoutReal neutral_energy_flow_to_core = 2.0 * nncore * tncore * v_th * dacore; // Incident energy
            BoutReal ionise_energy_flow = neutral_energy_flow_to_core * multiplier;

            // BoutReal core_neutral_energy_sink_flow = 0.0;
            // if (radial_energy_outflow(mesh->xstart-1, Viy, iz) > 0.0) {
            //   core_neutral_energy_sink_flow = 
            //     multiplier * radial_energy_outflow(mesh->xstart, iy, iz);
            //   ionise_energy_flow = core_neutral_energy_sink_flow;
            // }

            // diagnostic
            // source must be in domain
            channel.core_ion_density_source(mesh->xstart, iy, iz) = 
              ionise_particle_flow / volume;
            channel.core_ion_momentum_source(mesh->xstart, iy, iz) = 
              ionise_momentum_flow / volume;
            channel.core_ion_energy_source(mesh->xstart, iy, iz) = 
              ionise_energy_flow / volume;
            channel.core_neutral_density_sink(mesh->xstart, iy, iz) = 
              - ionise_particle_flow / volume;
            channel.core_neutral_momentum_sink(mesh->xstart, iy, iz) = 
              - ionise_momentum_flow / volume;
            channel.core_neutral_energy_sink(mesh->xstart, iy, iz) = 
              - ionise_energy_flow / volume;

            // save to solver
            ion_density_source(mesh->xstart, iy, iz) += 
              ionise_particle_flow / volume;
            ion_momentum_source(mesh->xstart, iy, iz) += 
              ionise_momentum_flow / volume;
            ion_energy_source(mesh->xstart, iy, iz) += 
              ionise_energy_flow / volume;
            neutral_density_source(mesh->xstart, iy, iz) -= 
              ionise_particle_flow / volume;
            neutral_momentum_source(mesh->xstart, iy, iz) -= 
              ionise_momentum_flow / volume;
            neutral_energy_source(mesh->xstart, iy, iz) -= 
              ionise_energy_flow / volume;
          }
        }
      }

      // Put the updated sources back into the state
      set<Field3D>(species_to["density_source"], ion_density_source);
      set<Field3D>(species_to["momentum_source"], ion_momentum_source);
      set<Field3D>(species_to["energy_source"], ion_energy_source);
      set<Field3D>(species_from["density_source"], neutral_density_source);
      set<Field3D>(species_from["momentum_source"], neutral_momentum_source);
      set<Field3D>(species_from["energy_source"], neutral_energy_source);
    }
  }
}

void NeutralBoundary::outputVars(Options& state) {

  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Cs0 = get<BoutReal>(state["Cs0"]);
  BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

  if (diagnose) {

    // Save particle and energy source for the species created during recycling

    // Target recycling

    if ((sol) or (pfr)) {
      set_with_attrs(
          state[{std::string("E") + name + std::string("_wall_refl")}],
          wall_energy_source,
          {{"time_dimension", "t"},
           {"units", "W m^-3"},
           {"conversion", Pnorm * Omega_ci},
           {"standard_name", "energy source"},
           {"long_name", std::string("Wall reflection energy source of ") + name},
           {"source", "neutral_boundary"}});
    }

    set_with_attrs(
        state[{std::string("E") + name + std::string("_target_refl")}],
        target_energy_source,
        {{"time_dimension", "t"},
         {"units", "W m^-3"},
         {"conversion", Pnorm * Omega_ci},
         {"standard_name", "energy source"},
         {"long_name", std::string("Wall reflection energy source of ") + name},
         {"source", "neutral_boundary"}});

    if (core_ionising){
      // Save particle and energy sink for neutrals
      for (const auto& channel : channels) {
        set_with_attrs(
          state[{std::string("N") + channel.from + std::string("_core_sink")}],
          channel.core_neutral_density_sink,
          {{"time_dimension", "t"},
           {"units", "m^-3 s^-1"},
           {"conversion", Nnorm * Omega_ci},
           {"standard_name", "particle sink"},
           {"long_name", std::string("Core neutral particle sink of ") + channel.from},
           {"source", "neutral_boundary"}});
        set_with_attrs(
          state[{std::string("NV") + channel.from + std::string("_core_sink")}],
          channel.core_neutral_momentum_sink,
          {{"time_dimension", "t"},
           {"units", "kg m^-2 s^-2"},
           {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
           {"standard_name", "momentum sink"},
           {"long_name", std::string("Core neutral momentum sink of ") + channel.from},
           {"source", "neutral_boundary"}});
        set_with_attrs(
          state[{std::string("E") + channel.from + std::string("_core_sink")}],
            channel.core_neutral_energy_sink,
            {{"time_dimension", "t"},
             {"units", "W m^-3"},
             {"conversion", Pnorm * Omega_ci},
            {"standard_name", "energy sink"},
            {"long_name", std::string("Core neutral energy sink of ") + channel.from},
            {"source", "neutral_boundary"}});
        set_with_attrs(
            state[{std::string("N") + channel.to + std::string("_core_source")}],
            channel.core_ion_density_source,
            {{"time_dimension", "t"},
            {"units", "m^-3 s^-1"},
            {"conversion", Nnorm * Omega_ci},
            {"standard_name", "particle source"},
            {"long_name", std::string("Core ionising particle source of ") + channel.to},
            {"source", "neutral_boundary"}});
        set_with_attrs(
            state[{std::string("NV") + channel.to + std::string("_core_source")}],
            channel.core_ion_momentum_source,
            {{"time_dimension", "t"},
            {"units", "kg m^-2 s^-2"},
            {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
            {"standard_name", "momentum source"},
            {"long_name", std::string("Core ionising momentum source of ") + channel.to},
            {"source", "neutral_boundary"}});
        set_with_attrs(
            state[{std::string("E") + channel.to + std::string("_core_source") }],
            channel.core_ion_energy_source,
            {{"time_dimension", "t"},
            {"units", "W m^-3"},
            {"conversion", Pnorm * Omega_ci},
            {"standard_name", "energy source"},
            {"long_name", std::string("Core ionising energy source of ") + channel.to},
            {"source", "neutral_boundary"}});
      }
    }
  }
}
