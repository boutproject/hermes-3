#pragma once
#ifndef AMJUEL_REACTION_H
#define AMJUEL_REACTION_H

#include <algorithm>
#include <filesystem>
#include <string>

#include "amjueldata.hxx"
#include "component.hxx"
#include "integrate.hxx"
#include "reaction.hxx"

/**
 * @brief 
 * 
 */
struct AmjuelReaction : public Reaction {
  AmjuelReaction(std::string name, std::string short_reaction_type,
                 std::string amjuel_lbl, std::string from_species, std::string to_species,
                 Options& alloptions)
      : Reaction(name, alloptions), amjuel_data(short_reaction_type, amjuel_lbl),
        amjuel_src(std::string("Amjuel ") + amjuel_lbl), from_species(from_species),
        short_reaction_type(short_reaction_type), to_species(to_species) {}

protected:
  // Store some strings for use in attribute docstrings
  const std::string amjuel_src;
  const std::string short_reaction_type;
  const std::string from_species;
  const std::string to_species;

  /// Functions to calculate Amjuel rates from underlying tables
  BoutReal eval_amjuel_fit(BoutReal T, BoutReal n,
                           const std::vector<std::vector<BoutReal>>& coeff_table);
  virtual BoutReal eval_sigma_v_E(BoutReal T, BoutReal n) override final;
  virtual BoutReal eval_sigma_v(BoutReal T, BoutReal n) override final;

  virtual void transform_additional(Options& state,
                                    Field3D& reaction_rate) override final;

private:
  const AmjuelData amjuel_data;
};

#endif // AMJUEL_REACTION_H
