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

    Field3D geometric_broadening{0.0};
    Field3D transport_broadening{0.0};
    Field3D flux_expansion{0.0};

    Field3D flux_tube_broadening{0.0};
};

namespace {
  RegisterComponent<FieldlineGeometry> register_fieldline_geometry("fieldline_geometry");
}

#endif // fieldline_geometry_H