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

    Field3D lambda_int{0.0};
    Field3D fieldline_radius{0.0};
    Field3D poloidal_magnetic_field{0.0};
    Field3D toroidal_magnetic_field{0.0};
    Field3D total_magnetic_field{0.0};
    Field3D pitch_angle{0.0};
    
    Field3D transport_broadening{0.0};
    Field3D flux_expansion{0.0};

    Field3D flux_tube_width{0.0};
    Field3D cell_poloidal_length{0.0};
    Field3D cell_side_area{0.0};
    Field3D cell_volume{0.0};

};

namespace {
  RegisterComponent<FieldlineGeometry> register_fieldline_geometry("fieldline_geometry");
}

#endif // fieldline_geometry_H