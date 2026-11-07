#pragma once

#include "component.hxx"
#include "generated_sources.hxx"

#include <bout/field_factory.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
#include <bout/output.hxx>
#include <bout/sys/expressionparser.hxx>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

class IonEnergyMmsCache {
public:
  explicit IonEnergyMmsCache(const Mesh& mesh, int variant)
      : variant(variant), ny(static_cast<size_t>(mesh.LocalNy)),
        nz(static_cast<size_t>(mesh.LocalNz)),
        cells(static_cast<size_t>(mesh.LocalNx) * ny * nz) {}

  const std::array<BoutReal, 7>& get(const bout::generator::Context& ctx) {
    ASSERT2(ctx.ix() >= 0 && ctx.jy() >= 0 && ctx.kz() >= 0);
    const auto index =
        (static_cast<size_t>(ctx.ix()) * ny + static_cast<size_t>(ctx.jy())) * nz
        + static_cast<size_t>(ctx.kz());
    ASSERT2(index < cells.size());
    auto& cell = cells[index];
    if (cell.time != ctx.t()) {
      cell.values =
          ion_energy_mms_test::sources(variant, ctx.x(), ctx.y(), ctx.z(), ctx.t());
      cell.time = ctx.t();
    }
    return cell.values;
  }

private:
  struct Cell {
    BoutReal time{std::numeric_limits<BoutReal>::quiet_NaN()};
    std::array<BoutReal, 7> values{};
  };
  int variant;
  size_t ny, nz;
  std::vector<Cell> cells;
};

/// Fast forcing for test builds only. BOUT owns MMS source addition,
/// initialisation and native errors; no physical species sources are set.
class IonEnergyMmsGenerator : public FieldGenerator {
public:
  IonEnergyMmsGenerator(size_t index, std::shared_ptr<IonEnergyMmsCache> cache)
      : index(index), cache(std::move(cache)) {}

  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> args) override {
    if (!args.empty()) {
      throw BoutException("Ion energy MMS generators take no arguments");
    }
    return std::make_shared<IonEnergyMmsGenerator>(index, cache);
  }

  double generate(const bout::generator::Context& ctx) override {
    return cache->get(ctx)[index];
  }

private:
  size_t index;
  std::shared_ptr<IonEnergyMmsCache> cache;
};

struct TestIonEnergyMms : public NamedComponent<TestIonEnergyMms> {
  TestIonEnergyMms(std::string name, Options& alloptions, Solver*)
      : NamedComponent(name, {}) {
    if (!alloptions["solver"]["mms"].as<bool>()) {
      throw BoutException("test_ion_energy_mms requires solver:mms=true");
    }
    const auto* mesh = bout::globals::mesh;
    if (mesh->periodicX) {
      throw BoutException("Ion-energy slab MMS requires nonperiodic X");
    }
    for (int x = mesh->xstart; x <= mesh->xend; ++x) {
      if (!mesh->periodicY(x)) {
        throw BoutException("Ion-energy slab MMS requires periodic Y at every X");
      }
    }
    const std::array<std::string, 7> variables{
        {"Ni", "Nj", "NVi", "NVj", "Pe", "Vort", "Eion"}};
    auto* factory = FieldFactory::get();
    const auto variant_name = alloptions["mms"]["variant"].as<std::string>();
    const int variant =
        variant_name == "skew" ? 1 : (variant_name == "viscosity" ? 2 : 0);
    if (variant == 0 && variant_name != "baseline") {
      throw BoutException("Unknown ion-energy MMS variant: {}", variant_name);
    }
    auto cache = std::make_shared<IonEnergyMmsCache>(*bout::globals::mesh, variant);
    // Cross-check BOUT expressions at interior, boundary and time-dependent
    // points on every MPI rank before selecting the compiled evaluation.
    for (size_t index = 0; index < variables.size(); ++index) {
      auto& section = alloptions[variables[index]];
      const auto expression = section["source"].as<std::string>();
      const auto analytic = factory->parse(expression, &section);
      for (const auto& point :
           std::array<std::array<BoutReal, 4>, 4>{{{{0.137, 1.83, 2.66, 0.73}},
                                                   {{0.71, 4.22, 0.91, 2.31}},
                                                   {{0.0, 2.41, 4.37, 1.17}},
                                                   {{1.0, 0.33, 5.71, 4.01}}}}) {
        const auto context = bout::generator::Context().set("x", point[0], "y", point[1],
                                                            "z", point[2], "t", point[3]);
        const auto reference = analytic->generate(context);
        const auto value = ion_energy_mms_test::sources(variant, point[0], point[1],
                                                        point[2], point[3])[index];
        if (!std::isfinite(value)
            || std::abs(reference - value) > 1e-12 * (1.0 + std::abs(reference))) {
          throw BoutException("Compiled/BOUT MMS source mismatch for {}: {} versus {}",
                              variables[index], value, reference);
        }
      }
      const std::string generator_name = "ion_energy_mms_" + variables[index];
      factory->addGenerator(generator_name,
                            std::make_shared<IonEnergyMmsGenerator>(index, cache));
      alloptions[name][variables[index] + "_analytic_source"].force(expression);
      // BOUT adds MMS sources after finally(), including evolved guard cells.
      // Vort's radial guards store phi midpoints, so physical vorticity forcing
      // must be restricted to the interior rather than evolving that wall state.
      const std::string source =
          index == 5 ? "where(x*(1-x)," + generator_name + ",0)" : generator_name;
      section["source"].force(source, "test_ion_energy_mms");
    }
    factory->cleanCache();
    output_info.write("Compiled ion-energy MMS sources agree with BOUT expressions\n");
  }

  static constexpr auto type = "test_ion_energy_mms";

private:
  void transform_impl(GuardedOptions&) override {}
};

namespace {
RegisterComponent<TestIonEnergyMms> register_test_ion_energy_mms;
}
