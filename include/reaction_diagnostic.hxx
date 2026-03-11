#pragma once
#ifndef REACTION_DIAGNOSTIC_H
#define REACTION_DIAGNOSTIC_H

#include <functional>
#include <string>

#include <bout/bout_enum_class.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>

#include "component.hxx"

BOUT_ENUM_CLASS(ReactionDiagnosticType, collision_freq, density_src, energy_src,
                energy_loss, momentum_src);

std::string toString(ReactionDiagnosticType diag_type);

static const std::map<ReactionDiagnosticType, std::string> state_labels = {
    {ReactionDiagnosticType::collision_freq, "collision_frequency"},
    {ReactionDiagnosticType::density_src, "density_source"},
    {ReactionDiagnosticType::energy_src, "energy_source"},
    {ReactionDiagnosticType::energy_loss, "energy_source"},
    {ReactionDiagnosticType::momentum_src, "momentum_source"}};

typedef std::function<Field3D(const Field3D&)> DiagnosticTransformerType;

// Some pre-defined functions to use as diagnostic transformers. Default is identity.
inline Field3D identity(const Field3D& src_fld) { return src_fld; };
inline Field3D negate(const Field3D& src_fld) { return -src_fld; };

struct ReactionDiagnostic {
  static inline std::string default_std_name(ReactionDiagnosticType type) {
    std::string standard_name;
    switch (type) {
    case ReactionDiagnosticType::collision_freq:
      standard_name = "collision frequency";
      break;
    case ReactionDiagnosticType::density_src:
      standard_name = "particle source";
      break;
    case ReactionDiagnosticType::energy_loss:
      standard_name = "radiation loss";
      break;
    case ReactionDiagnosticType::energy_src:
      standard_name = "energy transfer";
      break;
    case ReactionDiagnosticType::momentum_src:
      standard_name = "momentum transfer";
      break;
    default:
      throw BoutException(
          fmt::format("ReactionDiagnostic type [{:s}]] not set up!", toString(type)));
    }
    return standard_name;
  }

public:
  ReactionDiagnostic(const std::string& name, const std::string& long_name,
                     ReactionDiagnosticType type, const std::string& source,
                     DiagnosticTransformerType transformer = identity)
      : ReactionDiagnostic(name, long_name, type, source, default_std_name(type),
                           transformer){};

  ReactionDiagnostic(const std::string& name, const std::string& long_name,
                     ReactionDiagnosticType type, const std::string& source,
                     const std::string& standard_name,
                     DiagnosticTransformerType transformer)
      : long_name(long_name), name(name), source(source), standard_name(standard_name),
        transformer(transformer), type(type) {

    // Extract normalisation units from the root options
    Options& options = Options::root();
    const BoutReal Nnorm = get<BoutReal>(options["units"]["inv_meters_cubed"]);
    const BoutReal Tnorm = get<BoutReal>(options["units"]["eV"]);
    const BoutReal Omega_ci = 1 / get<BoutReal>(options["units"]["seconds"]);

    switch (type) {
    case ReactionDiagnosticType::collision_freq: {
      this->conversion = Omega_ci;
      this->units = "s^-1";
      break;
    }
    case ReactionDiagnosticType::density_src: {
      this->conversion = Nnorm * Omega_ci;
      this->units = "m^-3 s^-1";
      break;
    }
    case ReactionDiagnosticType::energy_loss: // fall through
    case ReactionDiagnosticType::energy_src: {
      const BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation
      this->conversion = Pnorm * Omega_ci;
      this->units = "W / m^3";
      break;
    }
    case ReactionDiagnosticType::momentum_src: {
      const BoutReal Cs0 = std::sqrt(SI::qe * Tnorm / SI::Mp);
      this->conversion = SI::Mp * Nnorm * Cs0 * Omega_ci;
      this->units = "kg m^-2 s^-2";
      break;
    }
    default:
      throw BoutException("ReactionDiagnostic type not set up!");
    }
    this->data.allocate();
  }

  // No copies
  ReactionDiagnostic(const ReactionDiagnostic&) = delete;
  // No copy-assignment
  ReactionDiagnostic& operator=(const ReactionDiagnostic&) = delete;
  // No move-assignment
  ReactionDiagnostic& operator=(ReactionDiagnostic&&) = delete;
  // Allow moves (needed for std::map::insert)
  ReactionDiagnostic(ReactionDiagnostic&&) = default;
  ~ReactionDiagnostic() = default;

  void add_to_state(Options& state) {
    set_with_attrs(state[this->name], this->data,
                   {{"time_dimension", "t"},
                    {"units", this->units},
                    {"conversion", this->conversion},
                    {"standard_name", this->standard_name},
                    {"long_name", this->long_name},
                    {"source", this->source}});
  }

  std::string get_name() const { return this->name; }

  void set_data(const Field3D& data) { this->data = data; }

  Field3D transform(Field3D src_fld) { return this->transformer(src_fld); }

private:
  BoutReal conversion;
  Field3D data;
  const std::string long_name;
  const std::string name;
  const std::string source;
  const std::string standard_name;
  const DiagnosticTransformerType transformer;
  const ReactionDiagnosticType type;
  std::string units;
};

#endif
