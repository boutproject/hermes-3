#ifndef TEST_AMJUEL_REACTIONS_H__
#define TEST_AMJUEL_REACTIONS_H__

#include "test_reactions.hxx"

#include "../../include/amjuel_hyd_ionisation.hxx"
#include "../../include/amjuel_hyd_recombination.hxx"

/**
 * @brief Base fixture for tests of AmjuelReaction subclasses.
 *
 * @tparam RTYPE The reaction class; must derive from AmjuelReaction.
 */
template <typename RTYPE>
class AmjuelReactionTest : public ReactionTest<RTYPE> {

  static_assert(std::is_base_of<AmjuelReaction, RTYPE>(),
                "Template arg to AmjuelReactionTest must derive from AmjuelReaction");

public:
  AmjuelReactionTest(std::string lbl, std::string reaction_str, std::string isotope)
      : ReactionTest<RTYPE>(lbl, reaction_str), isotope(isotope) {}

protected:
  std::string isotope;

  Options generate_state() override {
    std::string ion = isotope + "+";

    const std::map<std::string, BoutReal> sp_masses = {
        {"h", 1.0}, {"d", 2.0}, {"t", 3.0}, {"e", 1. / 1836}};
    const BoutReal atom_charge = 0.0;
    const BoutReal ion_charge = 1.0;

    Options state{{"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"species",
                   {{"e",
                     {{"AA", sp_masses.at("e")},
                      {"density", 1.0},
                      {"temperature", 1.0},
                      {"velocity", 1.0},
                      {"density_source", 0.0},
                      {"momentum_source", 0.0},
                      {"energy_source", 0.0}}},
                    {isotope,
                     {{"AA", sp_masses.at(isotope)},
                      {"charge", atom_charge},
                      {"density", 1.0},
                      {"temperature", 1.0},
                      {"velocity", 1.0},
                      {"density_source", 0.0},
                      {"momentum_source", 0.0},
                      {"energy_source", 0.0}}},
                    {ion,
                     {{"AA", sp_masses.at(isotope)},
                      {"charge", ion_charge},
                      {"density", 1.0},
                      {"temperature", 1.0},
                      {"velocity", 1.0},
                      {"density_source", 0.0},
                      {"momentum_source", 0.0},
                      {"energy_source", 0.0}}}}}};

    // Overwrite n_e, n_T, n_ion with functions that vary linearly along one axis

    // Density and Temperature ranges (log vals)
    constexpr BoutReal logn_min = 14, logn_max = 22;
    constexpr BoutReal logT_min = -1, logT_max = 4;

    state["species"][ion]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logn_min, logn_max, "x", this->nx), &state, mesh);
    state["species"]["e"]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logn_min, logn_max, "y", this->ny), &state, mesh);
    state["species"]["e"]["temperature"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logT_min, logT_max, "z", this->nz), &state, mesh);

    return state;
  }
};

class AmjuelHRecTest : public AmjuelReactionTest<AmjuelHydRecombinationIsotope<'h'>> {
public:
  AmjuelHRecTest()
      : AmjuelReactionTest<AmjuelHydRecombinationIsotope<'h'>>("Hrec", "h+ + e -> h",
                                                               "h") {}
};

class AmjuelDRecTest : public AmjuelReactionTest<AmjuelHydRecombinationIsotope<'d'>> {
public:
  AmjuelDRecTest()
      : AmjuelReactionTest<AmjuelHydRecombinationIsotope<'d'>>("Drec", "d+ + e -> d",
                                                               "d") {}
};

class AmjuelTRecTest : public AmjuelReactionTest<AmjuelHydRecombinationIsotope<'t'>> {
public:
  AmjuelTRecTest()
      : AmjuelReactionTest<AmjuelHydRecombinationIsotope<'t'>>("Trec", "t+ + e -> t",
                                                               "t") {}
};

class AmjuelHIznTest : public AmjuelReactionTest<AmjuelHydIonisationIsotope<'h'>> {
public:
  AmjuelHIznTest()
      : AmjuelReactionTest<AmjuelHydIonisationIsotope<'h'>>("Hizn", "h + e -> h+ + 2e",
                                                            "h") {}
};

class AmjuelDIznTest : public AmjuelReactionTest<AmjuelHydIonisationIsotope<'d'>> {
public:
  AmjuelDIznTest()
      : AmjuelReactionTest<AmjuelHydIonisationIsotope<'d'>>("Dizn", "d + e -> d+ + 2e",
                                                            "d") {}
};

class AmjuelTIznTest : public AmjuelReactionTest<AmjuelHydIonisationIsotope<'t'>> {
public:
  AmjuelTIznTest()
      : AmjuelReactionTest<AmjuelHydIonisationIsotope<'t'>>("Tizn", "t + e -> t+ + 2e",
                                                            "t") {}
};

#endif