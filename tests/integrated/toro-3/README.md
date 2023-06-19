Toro test problem #3
====================

Solves 1D adiabatic fluid equations, evolving the density
(n), pressure (p) and parallel momentum (nv).

Continuity (density) equation:

    dn/dt + Div(n v) = 0

Pressure equation:

    dp/dt + Div(p v) = - (gamma-1) * p * Div(v)

Momentum equation:

    d(nv)/dt + Div(nv v) = -Grad(p)

This simulation contains a strong shock close to a contact
discontinuity.

Left state: n = 1.0, v = 0, p = 1000

Right state: n = 1.0, v = 0, p = 0.01

The simulation runs to time t = 0.04
