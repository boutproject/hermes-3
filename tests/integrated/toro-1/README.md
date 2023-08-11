Toro test problem #1
====================

Solves 1D adiabatic fluid equations, evolving the density
(n), pressure (p) and parallel momentum (nv).

Continuity (density) equation:

    dn/dt + Div(n v) = 0

Pressure equation:

    dp/dt + Div(p v) = - (gamma-1) * p * Div(v)

Momentum equation:

    d(nv)/dt + Div(nv v) = -Grad(p)

The initial conditions are similar to the Sod shock problem,
but the left state is moving to the right.

Left state: n = 1.0, v = 0.75, p = 1.0

Right state: n = 0.125, v = 0, p = 0.1

The simulation runs to time t = 0.8
