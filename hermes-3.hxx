/*
    Copyright Hermes-3 contributors, 2016-2025
              email: dudson2@llnl.gov

    This file is part of Hermes-3 (Hot ion, multifluid).

    Hermes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Hermes is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Hermes.  If not, see <http://www.gnu.org/licenses/>.
  
*/

class Hermes;

#ifndef HERMES_H
#define HERMES_H

#include <bout/physicsmodel.hxx>

#include "include/component_scheduler.hxx"

class Hermes : public PhysicsModel {
public:
  virtual ~Hermes() {}
protected:
  int init(bool restarting) override;
  int rhs(BoutReal t) override;
  int precon(BoutReal t, BoutReal gamma, BoutReal delta);

  /// Add variables to be written to the output file
  ///
  /// Adds units and then calls each component in turn
  void outputVars(Options& options) override;

  /// Add variables to restart file
  void restartVars(Options& options) override;
private:
  /// Organises and schedules model components
  std::unique_ptr<ComponentScheduler> scheduler;

  /// Stores the dimensional units
  Options units;

  /// The evolving state
  Options state;
  
  /// Input normalisation constants
  BoutReal Tnorm, Nnorm, Bnorm;
  /// Derived normalisation constants
  BoutReal Cs0, Omega_ci, rho_s0;
};


#endif // HERMES_H
