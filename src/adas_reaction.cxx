///
/// Read and interpolate OpenADAS rate coefficients
///
/// Parts of this code are based on atomic++ by Thomas Body (2017)
///    https://github.com/TBody/atomicpp
///    Copyright (c) 2017 Tom Body
/// That was based on the TBody/atomic1D code,
///    https://github.com/TBody/atomicpp
/// which is in turn based on the cfe316/atomic code
///   https://github.com/cfe316/atomic
///   Copyright (c) 2016 David Wagner <wagdav@gmail.com>
///   Copyright (c) 2016 Jacob Schwartz <jschwart@pppl.gov>
///

#include "../include/adas_reaction.hxx"
#include "../include/integrate.hxx"

#include "../external/json.hxx"

#include <fstream>
#include <iterator>

OpenADASRateCoefficient::OpenADASRateCoefficient(const std::string& filename, int level) {
  AUTO_TRACE();

  // Read the rate file
  std::ifstream json_file(filename);

  if (!json_file.good()) {
    throw BoutException("Could not read ADAS file '{}'", filename);
  }
  
  nlohmann::json data;
  json_file >> data;

  // Get the log coefficients
  std::vector<std::vector<std::vector<double>>> extract_log_coeff = data["log_coeff"];
  std::vector<double> extract_log_temperature = data["log_temperature"];
  std::vector<double> extract_log_density = data["log_density"];

  log_coeff = extract_log_coeff[level];
  log_temperature = extract_log_temperature;
  log_density = extract_log_density;

  // Store the range of parameters
  Tmin = pow(10, log_temperature.front());
  Tmax = pow(10, log_temperature.back());

  nmin = pow(10, log_density.front());
  nmax = pow(10, log_density.back());
}

namespace {
BoutReal clip(BoutReal value, BoutReal min, BoutReal max) {
  if (value < min)
    return min;
  if (value > max)
    return max;
  return value;
}

int get_high_index(const std::vector<BoutReal>& vec, BoutReal value) {
  ASSERT2(vec.size() > 1); // Need at least two elements

  // Iterator pointing to the first element greater than or equal to log10T
  auto high_it = lower_bound(vec.begin(), vec.end(), value);
  if (high_it == vec.end()) {
    --high_it; // Shift to the last element
  } else if (high_it == vec.begin()) {
    ++high_it; // Shift to the second element
  }
  int high_index = std::distance(vec.begin(), high_it);

  ASSERT2((high_index > 0) and (high_index < static_cast<int>(vec.size())));
  return high_index;
}

} // namespace

BoutReal OpenADASRateCoefficient::evaluate(BoutReal T, BoutReal n) {
  AUTO_TRACE();

  // Ensure that the inputs are in range
  BoutReal log10T = log10(clip(T, Tmin, Tmax));
  BoutReal log10n = log10(clip(n, nmin, nmax));

  // Get the upper index. Between 1 and size-1 inclusive
  int high_T_index = get_high_index(log_temperature, log10T);
  int high_n_index = get_high_index(log_density, log10n);

  // Construct the simple interpolation grid
  // Find weightings based on linear distance
  // w01 ------ w11    ne -> y
  //  | \     / |      |
  //  |  w(x,y) |    --/--Te -> x
  //  | /     \ |      |
  // w00 ------ w10

  int low_T_index = high_T_index - 1;

  BoutReal x = (log10T - log_temperature[low_T_index])
               / (log_temperature[high_T_index] - log_temperature[low_T_index]);

  int low_n_index = high_n_index - 1;

  BoutReal y = (log10n - log_density[low_n_index])
               / (log_density[high_n_index] - log_density[low_n_index]);

  BoutReal eval_log_coef = (log_coeff[low_T_index][low_n_index] * (1 - y)
                            + log_coeff[low_T_index][high_n_index] * y)
                               * (1 - x)
                           + (log_coeff[high_T_index][low_n_index] * (1 - y)
                              + log_coeff[high_T_index][high_n_index] * y)
                                 * x;
  return pow(10., eval_log_coef);
}

void OpenADAS::calculate_rates(Options& electron, Options& from_ion, Options& to_ion) {
  AUTO_TRACE();

  Field3D Ne = get<Field3D>(electron["density"]);
  Field3D Te = get<Field3D>(electron["temperature"]);

  Field3D N1 = get<Field3D>(from_ion["density"]);
  Field3D T1 = get<Field3D>(from_ion["temperature"]);
  Field3D V1 = get<Field3D>(from_ion["velocity"]);
  auto AA = get<BoutReal>(from_ion["AA"]);

  ASSERT1(AA == get<BoutReal>(to_ion["AA"]));

  Field3D reaction_rate = cellAverage(
      [&](BoutReal ne, BoutReal n1, BoutReal te) {
        return ne * n1 * rate_coef.evaluate(te * Tnorm, ne * Nnorm) * Nnorm / FreqNorm;
      },
      Ne.getRegion("RGN_NOBNDRY"))(Ne, N1, Te);

  // Particles
  subtract(from_ion["density_source"], reaction_rate);
  add(to_ion["density_source"], reaction_rate);

  // Momentum
  Field3D momentum_exchange = reaction_rate * AA * V1;

  subtract(from_ion["momentum_source"], momentum_exchange);
  add(to_ion["momentum_source"], momentum_exchange);

  // Ion energy
  Field3D energy_exchange = reaction_rate * (3. / 2) * T1;
  subtract(from_ion["energy_source"], energy_exchange);
  add(to_ion["energy_source"], energy_exchange);

  // Electron energy loss (radiation, ionisation potential)
  Field3D energy_loss = cellAverage(
      [&](BoutReal ne, BoutReal n1, BoutReal te) {
        return ne * n1 * radiation_coef.evaluate(te * Tnorm, ne * Nnorm) * Nnorm
          / (Tnorm * FreqNorm);
      },
      Ne.getRegion("RGN_NOBNDRY"))(Ne, N1, Te);

  // Loss is reduced by heating
  energy_loss -= (electron_heating / Tnorm) * reaction_rate;

  subtract(electron["energy_source"], energy_loss);
}
