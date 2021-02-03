

#include <bout/fv_ops.hxx>
#include <derivs.hxx>
#include <difops.hxx>

#include "../include/neutral_mixed.hxx"

using bout::globals::mesh;

NeutralMixed::NeutralMixed(const std::string& name, Options& alloptions, Solver *solver) : name(name) {
  AUTO_TRACE();

  auto& options = alloptions[name];

  // Evolving variables e.g name is "h" or "h+"
  solver->add(Nn, std::string("N") + name);
  solver->add(Pn, std::string("P") + name);
  solver->add(NVn, std::string("NV") + name);

  sheath_ydown = options["sheath_ydown"]
                     .doc("Enable wall boundary conditions at ydown")
                     .withDefault<bool>(true);

  sheath_yup = options["sheath_yup"]
                   .doc("Enable wall boundary conditions at yup")
                   .withDefault<bool>(true);

  neutral_gamma =
      options["neutral_gamma"]
          .doc("Heat flux to the wall q = neutral_gamma * n * T * cs")
          .withDefault(5. / 4);

  nn_floor = options["nn_floor"]
                 .doc("A minimum density used when dividing NVn by Nn. "
                      "Normalised units.")
                 .withDefault(1e-4);

  precondition = options["precondition"]
                     .doc("Enable preconditioning in neutral model?")
                     .withDefault<bool>(true);
  
  if (precondition) {
    inv = std::unique_ptr<Laplacian>(Laplacian::create(&options["precon_laplace"]));

    inv->setInnerBoundaryFlags(INVERT_DC_GRAD | INVERT_AC_GRAD);
    inv->setOuterBoundaryFlags(INVERT_DC_GRAD | INVERT_AC_GRAD);

    inv->setCoefA(1.0);
  }

  // Optionally output time derivatives
  if (options["output_ddt"]
          .doc("Save derivatives to output?")
          .withDefault<bool>(false)) {
    bout::globals::dump.addRepeat(ddt(Nn), std::string("ddt(N") + name + std::string(")"));
    bout::globals::dump.addRepeat(ddt(Pn), std::string("ddt(P") + name + std::string(")"));
    bout::globals::dump.addRepeat(ddt(NVn), std::string("ddt(NV") + name + std::string(")"));
  }

  if (options["diagnose"].doc("Save additional diagnostics?").withDefault<bool>(false)) {
    bout::globals::dump.addRepeat(Sn, std::string("SN") + name);
    bout::globals::dump.addRepeat(Sp, std::string("SP") + name);
    bout::globals::dump.addRepeat(Snv, std::string("SNV") + name);
    Sn = Sp = Snv = 0.0;
  }

  AA = options["AA"].doc("Particle atomic mass. Proton = 1").withDefault(1.0);
}

void NeutralMixed::transform(Options &state) {
  AUTO_TRACE();

  mesh->communicate(Nn, Pn, NVn);

  Nn.clearParallelSlices();
  Pn.clearParallelSlices();
  NVn.clearParallelSlices();

  Nn = floor(Nn, 1e-8);
  Pn = floor(Pn, 1e-10);
  
  // Nnlim Used where division by neutral density is needed
  Nnlim = floor(Nn, nn_floor);
  Tn = Pn / Nn;
  // Tn = floor(Tn, 0.01 / Tnorm);
  Vn = NVn / (AA * Nn);
  Vnlim = NVn / (AA * Nnlim); // Neutral parallel velocity

  Pnlim = Nn * Tn;
  Pnlim.applyBoundary("neumann");

  /////////////////////////////////////////////////////
  // Boundary conditions
  TRACE("Neutral boundary conditions");

  if (sheath_ydown) {
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        // Free boundary (constant gradient) density
        BoutReal nnwall = 0.5 * (3. * Nn(r.ind, mesh->ystart, jz) -
                                 Nn(r.ind, mesh->ystart + 1, jz));
        if (nnwall < 0.0)
          nnwall = 0.0;

        BoutReal tnwall = Tn(r.ind, mesh->ystart, jz);

        Nn(r.ind, mesh->ystart - 1, jz) =
            2 * nnwall - Nn(r.ind, mesh->ystart, jz);

        // Zero gradient temperature, heat flux added later
        Tn(r.ind, mesh->ystart - 1, jz) = tnwall;

        // Set pressure consistent at the boundary
        // Pn(r.ind, mesh->ystart - 1, jz) =
        //     2. * nnwall * tnwall - Pn(r.ind, mesh->ystart, jz);

        // Zero-gradient pressure
        Pn(r.ind, mesh->ystart - 1, jz) = Pn(r.ind, mesh->ystart, jz);

        // No flow into wall
        Vn(r.ind, mesh->ystart - 1, jz) = -Vn(r.ind, mesh->ystart, jz);
        Vnlim(r.ind, mesh->ystart - 1, jz) = -Vnlim(r.ind, mesh->ystart, jz);
        NVn(r.ind, mesh->ystart - 1, jz) = -NVn(r.ind, mesh->ystart, jz);
      }
    }
  }

  if (sheath_yup) {
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        // Free boundary (constant gradient) density
        BoutReal nnwall = 0.5 * (3. * Nn(r.ind, mesh->yend, jz) -
                                 Nn(r.ind, mesh->yend - 1, jz));
        if (nnwall < 0.0)
          nnwall = 0.0;

        BoutReal tnwall = Tn(r.ind, mesh->yend, jz);

        Nn(r.ind, mesh->yend + 1, jz) = 2 * nnwall - Nn(r.ind, mesh->yend, jz);

        // Zero gradient temperature, heat flux added later
        Tn(r.ind, mesh->yend + 1, jz) = tnwall;

        // Set pressure consistent at the boundary
        // Pn(r.ind, mesh->yend + 1, jz) =
        //     2. * nnwall * tnwall - Pn(r.ind, mesh->yend, jz);

        // Zero-gradient pressure
        Pn(r.ind, mesh->yend + 1, jz) = Pn(r.ind, mesh->yend, jz);

        // No flow into wall
        Vn(r.ind, mesh->yend + 1, jz) = -Vn(r.ind, mesh->yend, jz);
        Vnlim(r.ind, mesh->yend + 1, jz) = -Vnlim(r.ind, mesh->yend, jz);
        NVn(r.ind, mesh->yend + 1, jz) = -NVn(r.ind, mesh->yend, jz);
      }
    }
  }
  

  // Set values in the state
  auto& localstate = state["species"][name];
  set(localstate["density"], Nn);
  set(localstate["AA"], AA); // Atomic mass
  set(localstate["pressure"], Pn);
  set(localstate["momentum"], NVn);
  set(localstate["velocity"], Vn);
  set(localstate["temperature"], Tn);
}

void NeutralMixed::finally(const Options &state) {
  AUTO_TRACE();
  auto& localstate = state["species"][name];

  ///////////////////////////////////////////////////////
  // Calculate cross-field diffusion from collision frequency
  //
  //
  BoutReal neutral_lmax =
      0.1 / get<BoutReal>(state["units"]["meters"]); // Normalised length
  Field3D Rnn = Nn * sqrt(Tn / AA) / neutral_lmax; // Neutral-neutral collisions
  
  if (localstate.isSet("collision_frequency")) {
    Dnn = Pnlim / (get<Field3D>(localstate["collision_frequency"]) + Rnn);
  } else {
    Dnn = Pnlim / Rnn;
  }

  mesh->communicate(Dnn);
  Dnn.clearParallelSlices();
  Dnn.applyBoundary("dirichlet_o2");

  // Apply a Dirichlet boundary condition to all the coefficients
  // used in diffusion operators. This is to ensure that the flux through
  // the boundary is zero.
  Field3D DnnPn = Dnn * Pn;
  DnnPn.applyBoundary("dirichlet_o2");
  Field3D DnnNn = Dnn * Nn;
  DnnNn.applyBoundary("dirichlet_o2");
  Field3D DnnNVn = Dnn * NVn;
  DnnNVn.applyBoundary("dirichlet_o2");

  if (sheath_ydown) {
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        Dnn(r.ind, mesh->ystart - 1, jz) = -Dnn(r.ind, mesh->ystart, jz);
        DnnPn(r.ind, mesh->ystart - 1, jz) = -DnnPn(r.ind, mesh->ystart, jz);
        DnnNn(r.ind, mesh->ystart - 1, jz) = -DnnNn(r.ind, mesh->ystart, jz);
        DnnNVn(r.ind, mesh->ystart - 1, jz) = -DnnNVn(r.ind, mesh->ystart, jz);
      }
    }
  }

  if (sheath_yup) {
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        Dnn(r.ind, mesh->yend + 1, jz) = -Dnn(r.ind, mesh->yend, jz);
        DnnPn(r.ind, mesh->yend + 1, jz) = -DnnPn(r.ind, mesh->yend, jz);
        DnnNn(r.ind, mesh->yend + 1, jz) = -DnnNn(r.ind, mesh->yend, jz);
        DnnNVn(r.ind, mesh->yend + 1, jz) = -DnnNVn(r.ind, mesh->yend, jz);
      }
    }
  }

  // Logarithms used to calculate perpendicular velocity
  // V_perp = -Dnn * ( Grad_perp(Nn)/Nn + Grad_perp(Tn)/Tn )
  //
  // Grad(Pn) / Pn = Grad(Tn)/Tn + Grad(Nn)/Nn
  //               = Grad(logTn + logNn)
  // Field3D logNn = log(Nn);
  // Field3D logTn = log(Tn);

  Field3D logPnlim = log(Pnlim);
  logPnlim.applyBoundary("neumann");

  // Sound speed appearing in Lax flux for advection terms
  Field3D sound_speed = sqrt(Tn * (5. / 3) / AA);
  
  /////////////////////////////////////////////////////
  // Neutral density
  TRACE("Neutral density");
  ddt(Nn) = -FV::Div_par(Nn, Vnlim, sound_speed) // Advection
            + FV::Div_a_Laplace_perp(DnnNn, logPnlim) // Perpendicular diffusion
    ;

  if (localstate.isSet("density_source")) {
    Sn = get<Field3D>(localstate["density_source"]);
    ddt(Nn) += Sn;
  }
  
  /////////////////////////////////////////////////////
  // Neutral momentum
  TRACE("Neutral momentum");

  // Relationship between heat conduction and viscosity for neutral
  // gas Chapman, Cowling "The Mathematical Theory of Non-Uniform
  // Gases", CUP 1952 Ferziger, Kaper "Mathematical Theory of
  // Transport Processes in Gases", 1972
  // eta_n = (2. / 5) * kappa_n;

  ddt(NVn) =
      -FV::Div_par(NVn, Vnlim, sound_speed)          // Momentum flow
      - Grad_par(Pnlim)                              // Pressure gradient
      + FV::Div_a_Laplace_perp(DnnNVn, logPnlim)     // Perpendicular diffusion
      + FV::Div_a_Laplace_perp((2. / 5) * DnnNn, Vn) // Perpendicular viscosity
      + FV::Div_par_K_Grad_par((2. / 5) * DnnNn, Vn) // Parallel viscosity
      ;

  if (localstate.isSet("momentum_source")) {
    Snv = get<Field3D>(localstate["momentum_source"]);
    ddt(NVn) += Snv;
  }

  /////////////////////////////////////////////////////
  // Neutral pressure
  TRACE("Neutral pressure");

  ddt(Pn) = -FV::Div_par(Pn, Vnlim, sound_speed) // Advection
            - (2. / 3) * Pnlim * Div_par(Vn)     // Compression
            + FV::Div_a_Laplace_perp(DnnPn, logPnlim) // Perpendicular diffusion
            + FV::Div_a_Laplace_perp(DnnNn, Tn)       // Conduction
            + FV::Div_par_K_Grad_par(DnnNn, Tn)       // Parallel conduction
      ;

  if (localstate.isSet("energy_source")) {
    Sp = (2. / 3) * get<Field3D>(localstate["energy_source"]);
    ddt(Pn) += Sp;
  }
  
  BOUT_FOR(i, Pn.getRegion("RGN_ALL")) {
    if ((Pn[i] < 1e-9) && (ddt(Pn)[i] < 0.0)) {
      ddt(Pn)[i] = 0.0;
    }
    if ((Nn[i] < 1e-7) && (ddt(Nn)[i] < 0.0)) {
      ddt(Nn)[i] = 0.0;
    }
  }
}

void NeutralMixed::annotate(Options &state) {
  auto& localstate = state["species"][name];
  localstate["density"].attributes["time_dimension"] = "t";
  localstate["density"].attributes["normalisation"] = "inv_meters_cubed";
}

void NeutralMixed::precon(const Options &state, BoutReal gamma) {
  if (!precondition) {
    return;
  }

  // Neutral gas diffusion
  // Solve (1 - gamma*Dnn*Delp2)^{-1}

  inv->setCoefD(-gamma * Dnn);

  ddt(Nn) = inv->solve(ddt(Nn));
  ddt(NVn) = inv->solve(ddt(NVn));
  ddt(Pn) = inv->solve(ddt(Pn));
}
