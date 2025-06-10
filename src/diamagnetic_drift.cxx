#include <bout/fv_ops.hxx>
#include <bout/vecops.hxx>
#include <bout/yboundary_regions.hxx>

#include "../include/diamagnetic_drift.hxx"

using bout::globals::mesh;

DiamagneticDrift::DiamagneticDrift(std::string name, Options& alloptions,
                                   Solver* UNUSED(solver)) {

  // Get options for this component
  auto& options = alloptions[name];

  yboundary.init(options);

  bndry_flux =
      options["bndry_flux"].doc("Allow fluxes through boundary?").withDefault<bool>(true);

  diamag_form = options["diamag_form"]
    .doc("Form of diamagnetic drift: 0 = gradient; 1 = divergence")
    .withDefault(Coordinates::FieldMetric(1.0));


  // from https://iopscience.iop.org/article/10.1088/1361-6587/aaed7d/pdf
  // Appendix A8
  bracket_form = options["bracket_form"]
    .doc("Use form of the diamagnetic drift that uses the poisson brackets?")
    .withDefault<bool>(true);

  
  // Read curvature vector
  Curlb_B.covariant = false; // Contravariant
  if (mesh->get(Curlb_B, "bxcv")) {
    throw BoutException("Curvature vector not found in input");
  }

  Options& paralleltransform = Options::root()["mesh"]["paralleltransform"];
  if (paralleltransform.isSet("type") and
      paralleltransform["type"].as<std::string>() == "shifted") {
    Field2D I;
    if (mesh->get(I, "sinty")) {
      I = 0.0;
    }
    Curlb_B.z += I * Curlb_B.x;
  }

  // Normalise

  // Get the units
  const auto& units = alloptions["units"];
  BoutReal Bnorm = get<BoutReal>(units["Tesla"]);
  BoutReal Lnorm = get<BoutReal>(units["meters"]);

  if (mesh->isFci()) {
    // All coordinates (x,y,z) are dimensionless
    // -> e_x has dimensions of length
    Curlb_B.x *= SQ(Lnorm);
  } else {
    // Field-aligned (Clebsch) coordinates
    Curlb_B.x /= Bnorm;
  }
  Curlb_B.y *= SQ(Lnorm);
  Curlb_B.z *= SQ(Lnorm);

  Curlb_B *= 2. / mesh->getCoordinates()->Bxy;

  mesh->communicate(Curlb_B.y);

  // Set drift to zero through sheath boundaries.
  // Flux through those cell faces should be set by sheath.
  yboundary.iter_regions([&](auto& region) {
    for (auto& pnt : region) {
      pnt.ynext(Curlb_B.y) = -Curlb_B.y[pnt.ind()];
    }
  });

  //FCI specific
  
  if (mesh->isFci()) {
    const auto coord = mesh->getCoordinates();
    bracket_factor = sqrt(coord->g_22.withoutParallelSlices()) / (coord->J.withoutParallelSlices() * coord->Bxy);

    logB = log(coord->Bxy);
    logB.applyBoundary("neumann_o2");
    mesh->communicate(logB);
    logB.applyParallelBoundary("parallel_neumann_o2");
    
  } else {
    bracket_factor = 1.0;
  }

  
  
}

void DiamagneticDrift::transform(Options& state) {
  // Iterate through all subsections
  Options& allspecies = state["species"];

  for (auto& kv : allspecies.getChildren()) {
    Options& species = allspecies[kv.first]; // Note: Need non-const

    if (!(species.isSet("charge") and species.isSet("temperature")))
      continue; // Skip, go to next species

    // Calculate diamagnetic drift velocity for this species
    auto q = get<BoutReal>(species["charge"]);
    if (fabs(q) < 1e-5) {
      continue;
    }
    auto T = GET_VALUE(Field3D, species["temperature"]);

    // Diamagnetic drift velocity
    Vector3D vD = (T / q) * Curlb_B;

    if (IS_SET(species["density"])) {
      auto N = GET_VALUE(Field3D, species["density"]);
      
      if (bracket_form){
	Field3D jdia_bracket = 2 * bracket(logB, T*N/q, BRACKET_ARAKAWA) * bracket_factor;
	subtract(species["density_source"], jdia_bracket);
      } else {  
	// Divergence form: Div(n v_D)
	Field3D div_form = FV::Div_f_v(N, vD, bndry_flux);
	// Gradient form: Curlb_B dot Grad(N T / q)
	Field3D grad_form = Curlb_B * Grad(N * T / q);
	subtract(species["density_source"], diamag_form * div_form + (1. - diamag_form) * grad_form);
      }
      
    }

    if (IS_SET(species["pressure"])) {
      auto P = get<Field3D>(species["pressure"]);
      if (bracket_form){
	Field3D jdia_bracket = 2 * bracket(logB, T*P/q, BRACKET_ARAKAWA) * bracket_factor;
	subtract(species["energy_source"], (5. / 2) * jdia_bracket);
      } else {
	Field3D div_form = FV::Div_f_v(P, vD, bndry_flux);
	Field3D grad_form = Curlb_B * Grad(P * T / q);
	subtract(species["energy_source"], (5. / 2) * (diamag_form * div_form + (1. - diamag_form) * grad_form));
      }
    }

    if (IS_SET(species["momentum"])) {
      auto NV = get<Field3D>(species["momentum"]);

      if (bracket_form){
	Field3D jdia_bracket = 2 * bracket(logB, T*NV/q, BRACKET_ARAKAWA) * bracket_factor;
	subtract(species["momentum_source"], jdia_bracket);
      } else {		   
	Field3D div_form = FV::Div_f_v(NV, vD, bndry_flux);
	Field3D grad_form = Curlb_B * Grad(NV * T / q);
	subtract(species["momentum_source"], diamag_form * div_form + (1. - diamag_form) * grad_form);
      }
    }
  }
}
