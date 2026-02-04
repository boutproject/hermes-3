# 2D energy test

Evolves vorticity, electron density, electron and ion pressure in a 2D
(32 x 32) box with curvature. Includes kinematic viscosity and ion
polarisation compression term.

The simulation starts with a non-uniform pressure that is unstable to
interchange modes. Electron and ion thermal energy is converted to
perpendicular motion. Viscosity then dissipates the flow, converting
the flow energy mainly to ion thermal energy.

The expected result is that electron pressure at the end is lower than
the start, because it has been converted to perpendicular motion. The
ion pressure should be higher at the end than at the start, because
motion is converted to ion heat. The total (ion + electron) pressure
should be higher at the end than the start, because potential energy
in the initial unstable state has been converted to thermal energy.

This test checks the sign of the electron, ion and total pressure
changes.  It is very low resolution so does not check the values or
convergence.
