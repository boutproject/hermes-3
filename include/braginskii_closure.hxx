#pragma once
#ifndef BRAGINSKII_H
#define BRAGINSKII_H

#include <bout/field3d.hxx>

#include "component.hxx"

/// Meta-component to set up all components necessary for the
/// Braginskii closure: `braginskii_collisions`,
/// `braginskii_friction`, `braginskii_heat_exchange`,
/// `braginskii_conduction`, `braginskii_electron_viscosity`,
/// `braginskii_ion_viscosity`, and `braginskii_thermal_force`. Each
/// of these components will have the same name as its type.
class BraginskiiClosure : public Component {
public:
  /// @param alloptions Settings, which may include
  ///    - <name>
  ///      - electron_viscosity  : bool    Include electron viscosity? (default: true)
  ///      - ion_viscosity  : bool    Include ion viscosity? (default: true)
  ///      - thermal_force  : bool    Inlucde thermal force between species? (default: true)
  BraginskiiClosure(std::string name, Options& alloptions, Solver*) : Component({}) {
    Options& options = alloptions[name];
    electron_viscosity = options["electron_viscosity"]
                          .doc("Include electron viscosity terms?")
                          .withDefault<bool>(true);
    ion_viscosity = options["ion_viscosity"]
                          .doc("Include ion viscosity terms?")
                          .withDefault<bool>(true);
    thermal_force = options["thermal_force"]
                          .doc("Include thermal force terms?")
                          .withDefault<bool>(true);
  }

  virtual std::vector<ComponentInformation> additionalComponents() override {
    std::vector<ComponentInformation> result = {
        {"braginskii_collisions", "braginskii_collisions"},
        {"braginskii_friction", "braginskii_friction"},
        {"braginskii_heat_exchange", "braginskii_heat_exchange"},
        {"braginskii_conduction", "braginskii_conduction"}};
    if (electron_viscosity) {
      result.emplace_back("braginskii_electron_viscosity",
                          "braginskii_electron_viscosity");
    }
    if (ion_viscosity) {
      result.emplace_back("braginskii_ion_viscosity", "braginskii_ion_viscosity");
    }
    if (thermal_force) {
      result.emplace_back("braginskii_thermal_force", "braginskii_thermal_force");
    }
    return result;
  }

private:
  bool electron_viscosity; /// Whether to include electron viscosity terms
  bool ion_viscosity;      /// Whether to include ion viscosity terms
  bool thermal_force;      /// Whether to include thermal force terms

  /// Empty transform; all the work actually happens in the subcomponents.
  void transform_impl(GuardedOptions&) override {}
};

namespace {
RegisterComponent<BraginskiiClosure>
    registercomponentbraginskiiclosure("braginskii_closure");
}

#endif // BRAGINSKII_H
