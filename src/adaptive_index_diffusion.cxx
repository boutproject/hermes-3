#include "adaptive_index_diffusion.hxx"

#include <bout/constants.hxx>
#include <bout/mesh.hxx>
#include <bout/utils.hxx>

#include <algorithm>
#include <cmath>
#include <vector>

namespace {
std::vector<std::string> enabledSections(const Options& options) {
  if (!options.isSet("components") && options.getChildren().count("hermes") == 0) {
    return {};
  }
  const auto& schedule = options.isSet("components") ? options : options["hermes"];
  std::vector<std::string> sections;
  if (schedule.isSet("components")) {
    for (const auto& entry : strsplit(schedule["components"].as<std::string>(), ',')) {
      const auto name = trim(entry, " \t\r\n()");
      if (!name.empty()) {
        sections.push_back(name);
      }
    }
  }
  return sections;
}

std::vector<std::string> adaptiveSections(const Options& options) {
  std::vector<std::string> sections;
  for (const auto& name : enabledSections(options)) {
    const auto types = options.getChildren().count(name) && options[name].isSet("type")
                           ? options[name]["type"].as<std::string>()
                           : name;
    for (const auto& type : strsplit(types, ',')) {
      if (trim(type, " \t\r\n()") == AdaptiveIndexDiffusion::type) {
        sections.push_back(name);
      }
    }
  }
  return sections;
}

BoutReal activationWeight(BoutReal value, BoutReal threshold) {
  if (!std::isfinite(value)) {
    throw BoutException("Nonfinite density or temperature in adaptive index diffusion");
  }
  value = std::max(value, 0.0);
  if (value <= threshold) {
    const auto ratio = value / threshold;
    return 1.0 / (1.0 + ratio * ratio);
  }
  const auto ratio = threshold / value;
  return ratio * ratio / (1.0 + ratio * ratio);
}
} // namespace

bool usesAdaptiveIndexDiffusion(const Options& options) {
  for (const auto& name : adaptiveSections(options)) {
    if (options.getChildren().count(name) == 0) {
      continue;
    }
    const auto& section = options[name];
    const auto rx =
        section.isSet("rate_x") ? section["rate_x"].as<Field2D>() : Field2D{0.0};
    const auto rz =
        section.isSet("rate_z") ? section["rate_z"].as<Field2D>() : Field2D{0.0};
    if (max(rx, true, "RGN_NOBNDRY") > 0.0 || max(rz, true, "RGN_NOBNDRY") > 0.0) {
      return true;
    }
  }
  return false;
}

AdaptiveIndexDiffusion::AdaptiveIndexDiffusion(std::string name, Options& options,
                                               Solver*)
    : NamedComponent(name, {}) {
  if (adaptiveSections(options).size() > 1) {
    throw BoutException("Enable only one adaptive_index_diffusion component");
  }
  auto& opts = options[name];
  const auto& units = options["units"];
  frequency = 1.0 / units["seconds"].as<BoutReal>();
  norm_density = units["inv_meters_cubed"];
  norm_temperature = units["eV"];
  if (!(std::isfinite(frequency) && frequency > 0.0 && std::isfinite(norm_density)
        && norm_density > 0.0 && std::isfinite(norm_temperature)
        && norm_temperature > 0.0)) {
    throw BoutException(
        "Adaptive index diffusion requires positive finite normalisations");
  }
  rate_x = opts["rate_x"]
               .doc("Maximum X cell loss rate [s^-1], may vary in X and Y")
               .withDefault<Field2D>(Field2D{0.0})
           / frequency;
  rate_z = opts["rate_z"]
               .doc("Maximum Z cell loss rate [s^-1], may vary in X and Y")
               .withDefault<Field2D>(Field2D{0.0})
           / frequency;
  for (const auto& i : rate_x.getRegion("RGN_NOBNDRY")) {
    if (!(std::isfinite(rate_x[i]) && rate_x[i] >= 0.0 && std::isfinite(rate_z[i])
          && rate_z[i] >= 0.0)) {
      throw BoutException(
          "Adaptive index diffusion rates must be finite and nonnegative");
    }
  }
  density_trigger = opts.isSet("density_threshold");
  temperature_trigger = opts.isSet("temperature_threshold");
  if (density_trigger) {
    density_threshold =
        opts["density_threshold"].doc("Half-activation density [m^-3]").as<BoutReal>()
        / norm_density;
  }
  if (temperature_trigger) {
    temperature_threshold = opts["temperature_threshold"]
                                .doc("Half-activation ion temperature [eV]")
                                .as<BoutReal>()
                            / norm_temperature;
    temperature_species =
        opts["temperature_species"]
            .doc("Ion used for temperature activation; omit for one ion")
            .withDefault(std::string{});
  }
  if ((density_trigger && !(std::isfinite(density_threshold) && density_threshold > 0.0))
      || (temperature_trigger
          && !(std::isfinite(temperature_threshold) && temperature_threshold > 0.0))) {
    throw BoutException(
        "Adaptive index diffusion thresholds must be positive and finite");
  }
  enabled =
      max(rate_x, true, "RGN_NOBNDRY") > 0.0 || max(rate_z, true, "RGN_NOBNDRY") > 0.0;
  if (!density_trigger && !temperature_trigger && enabled) {
    throw BoutException("Supply a density or temperature threshold for index diffusion");
  }
  diagnose = opts["diagnose"]
                 .doc("Save adaptive diffusion rates and increments")
                 .withDefault<bool>(false);
  // Zero rates must not introduce dependencies, snapshots, communication or
  // equation increments. In particular, they must preserve the physical
  // components' scheduling and preconditioning order.
  if (!enabled) {
    return;
  }
  for (const auto& permission :
       {readIfSet("species:{all_species}:charge"),
        readOnly("species:{positive_ions}:density", Regions::Interior),
        readIfSet("species:{charged}:density", Regions::Interior),
        readIfSet("species:{charged}:AA"),
        readIfSet("species:{charged}:velocity", Regions::Interior),
        readIfSet("species:{charged}:index_diffusion_density", Regions::Interior),
        readIfSet("species:{charged}:index_diffusion_momentum", Regions::Interior),
        readIfSet("species:{charged}:index_diffusion_density_floor"),
        readWrite("species:{charged}:index_diffusion_density_rhs"),
        readWrite("species:{charged}:index_diffusion_momentum_rhs"),
        readWrite("species:{charged}:index_diffusion_kinetic_energy"),
        readWrite("species:{charged}:index_diffusion_kinetic_rhs"),
        readWrite("fields:index_diffusion_conductance_x"),
        readWrite("fields:index_diffusion_conductance_z"),
        readIfSet("fields:Apar", Regions::Interior),
        readIfSet("fields:Apar_flutter", Regions::Interior)}) {
    setPermissions(permission);
  }
  if (temperature_trigger) {
    setPermissions(readIfSet("species:{positive_ions}:temperature", Regions::Interior));
  }
  for (const auto& section_name : enabledSections(options)) {
    if (options.getChildren().count(section_name) == 0) {
      continue;
    }
    const auto& section = options[section_name];
    if (!section.isSet("charge") || section["charge"].as<BoutReal>() == 0.0) {
      continue;
    }
    for (const auto* key :
         {"low_n_diffuse_perp", "low_p_diffuse_perp", "low_T_diffuse_perp"}) {
      if (section.isSet(key) && section[key].as<bool>()) {
        throw BoutException("adaptive_index_diffusion cannot be combined with {}:{}",
                            section_name, key);
      }
    }
  }
}

void AdaptiveIndexDiffusion::transform_impl(GuardedOptions& state) {
  if (!enabled) {
    return;
  }
  balance.clear();
  density_rates.clear();
  momentum_rates.clear();
  kinetic_rates.clear();
  if (state["fields"].isSet("Apar") || state["fields"].isSet("Apar_flutter")) {
    throw BoutException(
        "Adaptive index diffusion currently supports electrostatic models");
  }
  auto allspecies = state["species"];
  std::vector<std::string> ions, charged;
  Field3D charge_density{0.0}, electron_solver_density{0.0};
  for (const auto& [name, sp] : allspecies.getChildren()) {
    if (!sp.isSet("charge") || get<BoutReal>(sp["charge"]) == 0.0) {
      continue;
    }
    const auto charge = get<BoutReal>(sp["charge"]);
    if (charge < 0.0 && name != "e") {
      throw BoutException("Adaptive index diffusion does not support negative ion {}",
                          name);
    }
    charged.push_back(name);
    if (charge > 0.0) {
      ions.push_back(name);
      const auto density = getNoBoundary<Field3D>(sp["density"]);
      if (!sp.isSet("index_diffusion_density")) {
        throw BoutException("Adaptive index diffusion requires evolve_density for ion {}",
                            name);
      }
      const auto raw = getNoBoundary<Field3D>(sp["index_diffusion_density"]);
      BOUT_FOR(i, charge_density.getRegion("RGN_NOBNDRY")) {
        charge_density[i] += charge * density[i];
        electron_solver_density[i] += charge * raw[i];
      }
    }
  }
  if (ions.empty()) {
    throw BoutException("Adaptive index diffusion requires a positive ion species");
  }
  Field3D temperature{0.0};
  if (temperature_trigger) {
    const auto selected = temperature_species.empty() && ions.size() == 1
                              ? ions.front()
                              : temperature_species;
    if (std::find(ions.begin(), ions.end(), selected) == ions.end()
        || !allspecies[selected].isSet("temperature")) {
      throw BoutException(
          "Select an available ion temperature_species for index diffusion");
    }
    temperature = getNoBoundary<Field3D>(allspecies[selected]["temperature"]);
  }
  activation = Field3D{0.0};
  BOUT_FOR(i, activation.getRegion("RGN_NOBNDRY")) {
    const auto wn =
        density_trigger ? activationWeight(charge_density[i], density_threshold) : 0.0;
    const auto wt = temperature_trigger
                        ? activationWeight(temperature[i], temperature_threshold)
                        : 0.0;
    activation[i] = wn + (1.0 - wn) * wt;
  }
  const auto faces = indexDiffusionFaces(rate_x, rate_z, activation);
  loss_rate = indexDiffusionLossRate(faces);
  set(state["fields"]["index_diffusion_conductance_x"], faces.x);
  set(state["fields"]["index_diffusion_conductance_z"], faces.z);

  // First form all actual solver-moment increments. Quasineutral electron
  // density inherits the charge-weighted ion increment with the same faces.
  Field3D electron_physical_density_rate = zeroFrom(charge_density);
  for (const auto& name : charged) {
    const auto sp = allspecies[name];
    if (sp.isSet("index_diffusion_density")) {
      density_rates[name] = applyIndexDiffusion(
          faces, getNoBoundary<Field3D>(sp["index_diffusion_density"]));
    } else if (name == "e") {
      Field3D raw = zeroFrom(electron_solver_density);
      const auto charge = get<BoutReal>(sp["charge"]);
      BOUT_FOR(i, raw.getRegion("RGN_NOBNDRY")) {
        raw[i] = -electron_solver_density[i] / charge;
      }
      density_rates[name] = applyIndexDiffusion(faces, raw);
    }
    if (sp.isSet("index_diffusion_momentum")) {
      momentum_rates[name] = applyIndexDiffusion(
          faces, getNoBoundary<Field3D>(sp["index_diffusion_momentum"]));
    } else {
      momentum_rates[name] = Field3D{0.0};
    }
  }
  for (const auto& name : ions) {
    const auto sp = allspecies[name];
    const auto raw = getNoBoundary<Field3D>(sp["index_diffusion_density"]);
    const auto charge = get<BoutReal>(sp["charge"]);
    BOUT_FOR(i, electron_physical_density_rate.getRegion("RGN_NOBNDRY")) {
      if (raw[i] > 0.0) {
        electron_physical_density_rate[i] += charge * density_rates.at(name)[i];
      }
    }
  }
  for (const auto& name : charged) {
    auto sp = allspecies[name];
    const auto density = getNoBoundary<Field3D>(sp["density"]);
    Field3D raw = sp.isSet("index_diffusion_density")
                      ? getNoBoundary<Field3D>(sp["index_diffusion_density"])
                      : zeroFrom(density);
    if (name == "e" && !sp.isSet("index_diffusion_density")) {
      const auto charge = get<BoutReal>(sp["charge"]);
      BOUT_FOR(i, raw.getRegion("RGN_NOBNDRY")) {
        raw[i] = -electron_solver_density[i] / charge;
      }
    }
    const auto& sn = density_rates.at(name);
    const auto& sm = momentum_rates.at(name);
    const bool quasineutral_electron =
        name == "e" && !sp.isSet("index_diffusion_density");
    const BoutReal electron_fraction =
        quasineutral_electron ? -1.0 / get<BoutReal>(sp["charge"]) : 0.0;
    Field3D kinetic = zeroFrom(density), kinetic_rate = zeroFrom(density);
    if (sp.isSet("velocity")) {
      const auto velocity = getNoBoundary<Field3D>(sp["velocity"]);
      const BoutReal mass = get<BoutReal>(sp["AA"]);
      if (!sp.isSet("index_diffusion_momentum")) {
        throw BoutException("Adaptive index diffusion requires evolve_momentum for {} "
                            "when velocity is set",
                            name);
      }
      const BoutReal density_floor = get<BoutReal>(sp["index_diffusion_density_floor"]);
      BOUT_FOR(i, kinetic.getRegion("RGN_NOBNDRY")) {
        const auto n = std::max(density[i], 0.0);
        const auto v = velocity[i];
        kinetic[i] = 0.5 * mass * n * v * v;
        // Hard density clipping has zero derivative on its negative branch.
        const auto physical_sn =
            quasineutral_electron ? electron_fraction * electron_physical_density_rate[i]
                                  : (raw[i] > 0.0 ? sn[i] : 0.0);
        const auto nlim = softFloor(n, density_floor);
        const auto floor_derivative = -std::expm1(-n / density_floor);
        kinetic_rate[i] = (n / nlim) * v * sm[i]
                          + 0.5 * mass * v * v * (1.0 - 2.0 * n * floor_derivative / nlim)
                                * physical_sn;
      }
    }
    kinetic_rates[name] = kinetic_rate;
    set(sp["index_diffusion_density_rhs"], sn);
    set(sp["index_diffusion_momentum_rhs"], sm);
    set(sp["index_diffusion_kinetic_energy"], kinetic);
    set(sp["index_diffusion_kinetic_rhs"], kinetic_rate);
    balance.set("particles_" + name, "index_diffusion", sn);
    balance.set("momentum_" + name, "index_diffusion", sm);
    balance.set("kinetic_" + name, "index_diffusion", kinetic_rate);
  }
}

void AdaptiveIndexDiffusion::outputVars(Options& state) {
  balance.output(state, "adaptive_index_diffusion", diagnose);
  if (!diagnose || !activation.isAllocated()) {
    return;
  }
  set_with_attrs(state["index_diffusion_activation"], activation,
                 {{"time_dimension", "t"},
                  {"units", "1"},
                  {"long_name", "Adaptive index diffusion activation"}});
  set_with_attrs(state["index_diffusion_loss_rate"], loss_rate,
                 {{"time_dimension", "t"},
                  {"units", "s^-1"},
                  {"conversion", frequency},
                  {"long_name", "Outgoing index diffusion cell rate"}});
  const auto output_rates = [&](const auto& rates, const std::string& prefix,
                                const std::string& units, BoutReal conversion) {
    for (const auto& [name, rate] : rates) {
      set_with_attrs(state[prefix + name], rate,
                     {{"time_dimension", "t"},
                      {"units", units},
                      {"conversion", conversion},
                      {"species", name},
                      {"source", "adaptive_index_diffusion"}});
    }
  };
  output_rates(density_rates, "index_diffusion_density_rhs_", "m^-3 s^-1",
               norm_density * frequency);
  const auto velocity_norm = std::sqrt(SI::qe * norm_temperature / SI::Mp);
  output_rates(momentum_rates, "index_diffusion_momentum_rhs_", "kg m^-2 s^-2",
               SI::Mp * norm_density * velocity_norm * frequency);
  output_rates(kinetic_rates, "index_diffusion_kinetic_rhs_", "W m^-3",
               SI::qe * norm_density * norm_temperature * frequency);
}
