#pragma once
#ifndef fieldline_geometry_H
#define fieldline_geometry_H

#include "component.hxx"
#include <bout/constants.hxx>

struct FieldlineGeometry : public Component {

  FieldlineGeometry(std::string, Options& options, Solver*);

  void transform(Options& state) override;
  void outputVars(Options& state) override;

  private:
    BoutReal Lnorm;
    bool diagnose;

    Field3D lpar{0.0};

    Field3D B_poloidal{0.0};
    Field3D B_toroidal{0.0};
    Field3D major_radius{0.0};
    Field3D effective_flux_exp{0.0};
    Field3D poloidal_SOL_width{0.0};
    BoutReal lambda_q_omp;

    Field3D B_total{0.0};
    Field3D pitch_angle{0.0};
    Field3D area_external{0.0};

};

namespace {
  RegisterComponent<FieldlineGeometry> register_fieldline_geometry("fieldline_geometry");
}

#endif // fieldline_geometry_H