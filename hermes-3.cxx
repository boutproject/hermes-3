/*

    Copyright Hermes-3 contributors, 2016-2025
              email: dudson2@llnl.gov

    This file is part of Hermes-3 (Hot ion, multifluid)

    Hermes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Hermes is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Hermes.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "hermes-3.hxx"
#include "revision.hxx"

#include "include/adas_carbon.hxx"
#include "include/adas_neon.hxx"
#include "include/adas_lithium.hxx"
#include "include/amjuel_helium.hxx"
#include "include/amjuel_hyd_ionisation.hxx"
#include "include/amjuel_hyd_recombination.hxx"
#include "include/anomalous_diffusion.hxx"
#include "include/classical_diffusion.hxx"
#include "include/binormal_stpm.hxx"
#include "include/collisions.hxx"
#include "include/diamagnetic_drift.hxx"
#include "include/electromagnetic.hxx"
#include "include/electron_force_balance.hxx"
#include "include/electron_viscosity.hxx"
#include "include/evolve_density.hxx"
#include "include/evolve_energy.hxx"
#include "include/evolve_momentum.hxx"
#include "include/evolve_pressure.hxx"
#include "include/fixed_density.hxx"
#include "include/fixed_fraction_ions.hxx"
#include "include/fixed_fraction_radiation.hxx"
#include "include/fixed_temperature.hxx"
#include "include/fixed_velocity.hxx"
#include "include/hydrogen_charge_exchange.hxx"
#include "include/ion_viscosity.hxx"
#include "include/ionisation.hxx"
#include "include/isothermal.hxx"
#include "include/neutral_boundary.hxx"
#include "include/neutral_mixed.hxx"
#include "include/neutral_parallel_diffusion.hxx"
#include "include/noflow_boundary.hxx"
#include "include/polarisation_drift.hxx"
#include "include/quasineutral.hxx"
#include "include/recycling.hxx"
#include "include/relax_potential.hxx"
#include "include/scale_timederivs.hxx"
#include "include/set_temperature.hxx"
#include "include/sheath_boundary.hxx"
#include "include/sheath_boundary_insulating.hxx"
#include "include/sheath_boundary_simple.hxx"
#include "include/sheath_closure.hxx"
#include "include/simple_conduction.hxx"
#include "include/snb_conduction.hxx"
#include "include/solkit_hydrogen_charge_exchange.hxx"
#include "include/solkit_neutral_parallel_diffusion.hxx"
#include "include/sound_speed.hxx"
#include "include/thermal_force.hxx"
#include "include/transform.hxx"
#include "include/upstream_density_feedback.hxx"
#include "include/temperature_feedback.hxx"
#include "include/detachment_controller.hxx"
#include "include/vorticity.hxx"
#include "include/zero_current.hxx"
#include "include/simple_pump.hxx"
#include <bout/constants.hxx>
#include <bout/boundary_factory.hxx>
#include <bout/boundary_op.hxx>
#include <bout/field_factory.hxx>

#include "include/recalculate_metric.hxx"

#if !BOUT_USE_METRIC_3D
// For standard 2D metrics,
// Hermes operators don't need parallel slices
BOUT_OVERRIDE_DEFAULT_OPTION("mesh:calcParallelSlices_on_communicate", false);
#endif

class DecayLengthBoundary : public BoundaryOp {
public:
  DecayLengthBoundary() : gen(nullptr) {}
  DecayLengthBoundary(BoundaryRegion* region,  std::shared_ptr<FieldGenerator> g)
    : BoundaryOp(region), gen(std::move(g)) {}

  using BoundaryOp::clone;
  /// Create a copy of this boundary condition
  /// This is called by the Boundary Factory
  BoundaryOp* clone(BoundaryRegion* region, const std::list<std::string>& args) override {
    std::shared_ptr<FieldGenerator> newgen;
    if (!args.empty()) {
      // First argument should be an expression
      newgen = FieldFactory::get()->parse(args.front());
    }
    return new DecayLengthBoundary(region, newgen);
  }

  // Only implementing for Field3D, no time dependence
  void apply(Field3D& f) override {
    // Ensure that field and boundary are on the same mesh
    Mesh* mesh = bndry->localmesh;
    ASSERT1(mesh == f.getMesh());

    // Get cell radial length
    Coordinates *coord = mesh->getCoordinates();
    Field2D dx = coord->dx;
    Field2D g11 = coord->g11;
    Field2D dr = dx / sqrt(g11); // cell radial length. dr = dx/(Bpol * R) and g11 = (Bpol*R)**2 

    // Only implemented for cell centre quantities
    ASSERT1(f.getLocation() == CELL_CENTRE);

    // This loop goes over the first row of boundary cells (in X and Y)
    for (bndry->first(); !bndry->isDone(); bndry->next1d()) {
      for (int zk = 0; zk < mesh->LocalNz; zk++) { // Loop over Z points
        BoutReal decay_length = 3; // Default decay length is 3 normalised units (usually ~3mm)
        if (gen) {
          // Pick up the boundary condition setting from the input file
          // Must be specified in normalised units like the other BCs inputs
          decay_length = gen->generate(bout::generator::Context(bndry, zk, CELL_CENTRE, 0.0, mesh));
        }
        // Set value in inner guard cell f(bndry->x, bndry->y, zk)
        // using the final domain cell value f(bndry->x - bndry->bx, bndry->y - bndry->by, zk)
        // Note: (bx, by) is the direction into the boundary, so
        //    (1, 0)  X outer boundary (SOL)
        //    (-1, 0) X inner boundary (Core or PF)
        //    (0, 1)  Y upper boundary (outer lower target)
        //    (0, -1) Y lower boundary (inner lower target)
        
        // Distance between final cell centre and inner guard cell centre in normalised units
        BoutReal distance = 0.5 * (dr(bndry->x, bndry->y) + 
          dr(bndry->x - bndry->bx, bndry->y - bndry->by));

        // Exponential decay 
        f(bndry->x, bndry->y, zk) =
          f(bndry->x - bndry->bx, bndry->y - bndry->by, zk) * exp(-1 * distance / decay_length);

        // Set any remaining guard cells (i.e. the outer guards) to the same value
        // Should the outer guards have the decay continue, or just copy what the inners have?
        for (int i = 1; i < bndry->width; i++) {
          f(bndry->x + i * bndry->bx, bndry->y + i * bndry->by, zk) = f(bndry->x, bndry->y, zk);
        }
      }
    }
  }

  void apply(Field2D& f) override {
    throw BoutException("DecayLengthBoundary not implemented for Field2D");
  }
private:
  std::shared_ptr<FieldGenerator> gen; // Generator
};

int Hermes::init(bool restarting) {

  auto &options = Options::root()["hermes"];
  
  output.write("\nGit Version of Hermes: {:s}\n", hermes::version::revision);
  options["revision"] = hermes::version::revision;
  options["revision"].setConditionallyUsed();

  output.write("Slope limiter: {}\n", hermes::limiter_typename);
  options["slope_limiter"] = hermes::limiter_typename;
  options["slope_limiter"].setConditionallyUsed();

  // Choose normalisations
  Tnorm = options["Tnorm"].doc("Reference temperature [eV]").withDefault(100.);
  Nnorm = options["Nnorm"].doc("Reference density [m^-3]").withDefault(1e19);
  Bnorm = options["Bnorm"].doc("Reference magnetic field [T]").withDefault(1.0);

  Cs0 = sqrt(SI::qe * Tnorm / SI::Mp); // Reference sound speed [m/s]
  Omega_ci = SI::qe * Bnorm / SI::Mp;  // Ion cyclotron frequency [1/s]
  rho_s0 = Cs0 / Omega_ci;             // Length scale [m]

  // Put normalisation quantities into an Options to use later
  units["inv_meters_cubed"] = Nnorm;
  units["eV"] = Tnorm;
  units["Tesla"] = Bnorm;
  units["seconds"] = 1./Omega_ci;
  units["meters"] = rho_s0;

  // Put into the options tree, so quantities can be normalised
  // when creating components
  Options::root()["units"] = units.copy();
  Options::root()["units"].setConditionallyUsed();

  // Add the decay length boundary condition to the boundary factory
  // This will make it available as an input option
  // e.g. bndry_sol = decaylength(0.003 / rho_s0) sets up a decay length of 3mm
  BoundaryFactory::getInstance()->add(new DecayLengthBoundary(), "decaylength");

  /////////////////////////////////////////////////////////
  // Load metric tensor from the mesh, passing length and B
  // field normalisations
  TRACE("Loading metric tensor");

  if (options.isSet("loadmetric")) {
    throw BoutException("Error: The loadmetric option has been renamed to recalculate_metric.\n"
                        "Note: The default (true/false) is inverted for the new option.\n"
                        "Setting recalculate_metric=false (the default) uses the metric in the grid file.\n"
                        "Setting recalculate_metric=true loads Rxy, Bpxy etc from the grid file.\n"
                        "  This assumes an orthogonal coordinate system. See manual for details.");
  }

  if (options["recalculate_metric"]
          .doc("Load Rxy, Bpxy etc. to calculate an orthogonal metric?")
          .withDefault(false)) {
    recalculate_metric(rho_s0, Bnorm);

  } else {
    // Check that the grid file contains at least one metric tensor component
    // Note: Older grid files did not, so would silently default to the identity metric
    if (!(mesh->sourceHasVar("J") or
          mesh->sourceHasVar("g11") or
          mesh->sourceHasVar("g22") or
          mesh->sourceHasVar("g33") or
          mesh->sourceHasVar("g12") or
          mesh->sourceHasVar("g23") or
          mesh->sourceHasVar("g13") or
          mesh->sourceHasVar("g_11") or
          mesh->sourceHasVar("g_22") or
          mesh->sourceHasVar("g_33") or
          mesh->sourceHasVar("g_12") or
          mesh->sourceHasVar("g_23") or
          mesh->sourceHasVar("g_13"))) {
      throw BoutException("Grid input does not contain any metric components (J, g11, g_22 etc).\n"
                          "Set hermes:recalculate_metric=true to calculate from Rxy, Bpxy etc.\n"
                          "If the default identity metric is intended then set e.g. mesh:J=1\n");
    }

    if (options["normalise_metric"]
             .doc("Normalise input metric tensor? (assumes input is in SI units)")
             .withDefault<bool>(true)) {
      Coordinates *coord = mesh->getCoordinates();
      // To use non-orthogonal metric
      // Normalise
      coord->dx /= rho_s0 * rho_s0 * Bnorm;
      coord->Bxy /= Bnorm;
      // Metric is in grid file - just need to normalise
      coord->g11 /= SQ(Bnorm * rho_s0);
      coord->g22 *= SQ(rho_s0);
      coord->g33 *= SQ(rho_s0);
      coord->g12 /= Bnorm;
      coord->g13 /= Bnorm;
      coord->g23 *= SQ(rho_s0);

      coord->J *= Bnorm / rho_s0;

      coord->g_11 *= SQ(Bnorm * rho_s0);
      coord->g_22 /= SQ(rho_s0);
      coord->g_33 /= SQ(rho_s0);
      coord->g_12 *= Bnorm;
      coord->g_13 *= Bnorm;
      coord->g_23 /= SQ(rho_s0);

      coord->geometry(); // Calculate other metrics
    }
  }

  // Tell the components if they are restarting
  options["restarting"] = restarting;
  options["restarting"].setConditionallyUsed();

  TRACE("Creating components");

  // Create the components
  // Here options is passed as the scheduler configuration, so that
  // settings in [hermes] are used.
  // Options::root() is passed as the root of the component options, so that
  // individual components use their own sections, rather than subsections of [hermes].
  scheduler = ComponentScheduler::create(options, Options::root(), solver);

  // Preconditioner
  setPrecon((preconfunc)&Hermes::precon);

  return 0;
}

int Hermes::rhs(BoutReal time) {
  // Need to reset the state, since fields may be modified in transform steps
  state = Options();
  
  set(state["time"], time);
  state["units"] = units.copy();

  // Call all the components
  scheduler->transform(state);

  return 0;
}

/*!
 * Preconditioner. Solves the heat conduction
 *
 * @param[in] t  The simulation time
 * @param[in] gamma   Factor in front of the Jacobian in (I - gamma*J). Related
 * to timestep
 * @param[in] delta   Not used here
 */
int Hermes::precon(BoutReal t, BoutReal gamma, BoutReal UNUSED(delta)) {
  state["time"] = t;
  scheduler->precon(state, gamma);
  return 0;
}

void Hermes::outputVars(Options& options) {
  AUTO_TRACE();

  // Save the Hermes version in the output dump files
  options["HERMES_REVISION"].force(hermes::version::revision);

  options["HERMES_SLOPE_LIMITER"].force(hermes::limiter_typename);

  // Save normalisation quantities. These may be used by components
  // to calculate conversion factors to SI units

  set_with_attrs(options["Tnorm"], Tnorm, {
      {"units", "eV"},
      {"conversion", 1}, // Already in SI units
      {"standard_name", "temperature normalisation"},
      {"long_name", "temperature normalisation"}
    });
  set_with_attrs(options["Nnorm"], Nnorm, {
      {"units", "m^-3"},
      {"conversion", 1},
      {"standard_name", "density normalisation"},
      {"long_name", "Number density normalisation"}
    });
  set_with_attrs(options["Bnorm"], Bnorm, {
      {"units", "T"},
      {"conversion", 1},
      {"standard_name", "magnetic field normalisation"},
      {"long_name", "Magnetic field normalisation"}
    });
  set_with_attrs(options["Cs0"], Cs0, {
      {"units", "m/s"},
      {"conversion", 1},
      {"standard_name", "velocity normalisation"},
      {"long_name", "Sound speed normalisation"}
    });
  set_with_attrs(options["Omega_ci"], Omega_ci, {
      {"units", "s^-1"},
      {"conversion", 1},
      {"standard_name", "frequency normalisation"},
      {"long_name", "Cyclotron frequency normalisation"}
    });
  set_with_attrs(options["rho_s0"], rho_s0, {
      {"units", "m"},
      {"conversion", 1},
      {"standard_name", "length normalisation"},
      {"long_name", "Gyro-radius length normalisation"}
    });
  scheduler->outputVars(options);
}

void Hermes::restartVars(Options& options) {
  AUTO_TRACE();

  set_with_attrs(options["Tnorm"], Tnorm, {
      {"units", "eV"},
      {"conversion", 1}, // Already in SI units
      {"standard_name", "temperature normalisation"},
      {"long_name", "temperature normalisation"}
    });
  set_with_attrs(options["Nnorm"], Nnorm, {
      {"units", "m^-3"},
      {"conversion", 1},
      {"standard_name", "density normalisation"},
      {"long_name", "Number density normalisation"}
    });
  set_with_attrs(options["Bnorm"], Bnorm, {
      {"units", "T"},
      {"conversion", 1},
      {"standard_name", "magnetic field normalisation"},
      {"long_name", "Magnetic field normalisation"}
    });
  set_with_attrs(options["Cs0"], Cs0, {
      {"units", "m/s"},
      {"conversion", 1},
      {"standard_name", "velocity normalisation"},
      {"long_name", "Sound speed normalisation"}
    });
  set_with_attrs(options["Omega_ci"], Omega_ci, {
      {"units", "s^-1"},
      {"conversion", 1},
      {"standard_name", "frequency normalisation"},
      {"long_name", "Cyclotron frequency normalisation"}
    });
  set_with_attrs(options["rho_s0"], rho_s0, {
      {"units", "m"},
      {"conversion", 1},
      {"standard_name", "length normalisation"},
      {"long_name", "Gyro-radius length normalisation"}
    });
  scheduler->restartVars(options);
}

// Standard main() function
BOUTMAIN(Hermes);
