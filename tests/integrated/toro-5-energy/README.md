Toro test problem #5 (Energy form)
==================================

Solves 1D adiabatic fluid equations, evolving the density
(n), energy (E = 3/2p + 1/2 nv^2) and parallel momentum (nv).

Continuity (density) equation:

    dn/dt + Div(n v) = 0

Energy equation:

    dE/dt + Div([E + p] v) = 0

Momentum equation:

    d(nv)/dt + Div(nv v) = -Grad(p)

The initial conditions are the same as Toro's test problem #3
but moving to the left with a uniform speed so that the
contact discontinuity remains approximately stationary.

Left state: n = 1.0, v = -19.59745, p = 1000.0

Right state: n = 1.0, v = -19.59745, p = 0.01

The simulation runs to time t = 0.03
