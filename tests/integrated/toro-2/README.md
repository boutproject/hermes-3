Toro test problem #2
====================

Solves 1D adiabatic fluid equations, evolving the density
(n), pressure (p) and parallel momentum (nv).

Continuity (density) equation:

    dn/dt + Div(n v) = 0

Pressure equation:

    dp/dt + Div(p v) = - (gamma-1) * p * Div(v)

Momentum equation:

    d(nv)/dt + Div(nv v) = -Grad(p)

The initial condition has a diverging supersonic flow

Left state: n = 1.0, v = -2.0, p = 0.4

Right state: n = 1.0, v = 2.0, p = 0.4

The simulation runs to time t = 0.04
