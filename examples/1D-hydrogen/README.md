1D transport with hydrogen plasma
=================================

100% recycling of pure deuterium plasma. MAST-U like parameters,
similar to those used in SD1D simulations.

The solver settings are for the "beuler" backward Euler solver,
using the Jacobian coloring scheme. For comparison the "qn" subdirectory
contains an input for the Quasi-Newton method (matrix free, no
preconditioning).

After running both cases, running plot_convergence.py should
produce plots of the convergence of the density source against
RHS calls and simulation time, and the convergence of the RMS
time derivative against cumulative RHS calls.
