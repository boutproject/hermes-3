Toro test problem #4
====================

Solves 1D adiabatic fluid equations, evolving the density
(n), pressure (p) and parallel momentum (nv).

Continuity (density) equation:

    dn/dt + Div(n v) = 0

Pressure equation:

    dp/dt + Div(p v) = - (gamma-1) * p * Div(v)

Momentum equation:

    d(nv)/dt + Div(nv v) = -Grad(p)

This simulation develops two right-going shocks with a
contact discontinuity between them.

Left state: n = 5.99924, v = 19.5975, p = 460.894

Right state: n = 5.99242, v = -6.19633, p = 46.0950

The simulation runs to time t = 0.15
