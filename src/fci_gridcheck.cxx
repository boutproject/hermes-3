
#include "../include/fci_gridcheck.hxx"

#include <bout/constants.hxx>

#include <bout/mesh.hxx> 

using bout::globals::mesh;

FCIGridcheck::FCIGridcheck(std::string name, Options &alloptions, Solver *) {
  Options& options = alloptions[name];

  BoutReal Lnorm = alloptions["units"]["meters"]; // Length normalisation factor

  auto* coord = mesh->getCoordinates();

  auto B = coord->Bxy;

  auto g_22 = coord->g_22;

  auto J = coord->J;
  
  BOUT_FOR(i, B.getRegion("RGN_NOY")) {
    const auto iyp = i.yp();
    const auto iym = i.ym();

    BoutReal Bflux = B[i]*J[i] / sqrt(g_22[i]);
    BoutReal Bflux_up = B.yup()[iyp] * J.yup()[iyp] / sqrt(g_22.yup()[iyp]);
    BoutReal Bflux_down = B.ydown()[iym] * J.ydown()[iym] / sqrt(g_22.ydown()[iym]);

    BoutReal cor_up = Bflux / Bflux_up;
    BoutReal cor_down = Bflux / Bflux_down;

    J.yup()[iyp] *= cor_up;
    J.ydown()[iym] *= cor_down;
  }
  coord->J  = J;
  
}


void FCIGridcheck::transform(Options &state) {
  AUTO_TRACE();
}


void FCIGridcheck::outputVars(Options& state) {
  auto* coord = mesh->getCoordinates();

  set_with_attrs(state[std::string("J")], coord->J,
                     {{"units", "m^3"},
                      {"long_name", " Jacobian"},
                      {"source", "fci_gridcheck"}});
  
}
