.. _sec-tests:

Tests
=====

Sod shock
---------

``tests/integrates/sod-shock`` and ``tests/integrated/sod-shock-energy``

Euler equations in 1D. Starting from a state with a jump at the middle
of the domain.  Left state density, velocity and pressure are
:math:`\left(\rho_L, u_L, p_L\right) = \left(1.0, 0, 1.0\right)` Right
state :math:`\left(\rho_R, u_R, p_R\right) = \left(0.125, 0,
0.1\right)`. The result is shown in figure below at time :math:`t =
0.2` for different resolutions in a domain of length 1. The solid
black line is the analytic solution.

.. figure:: figs/sod_shock.png
   :name: sod_shock
   :alt:
   :width: 60%

When evolving pressure the position of the shock front lags the
analytic solution, with the pressure behind the front slightly too
high. This is a known consequence of solving the Euler equations in
non-conservative form. If instead we evolve energy (internal +
kinetic) then the result is much closer to the analytic solution.

.. figure:: figs/sod_shock_energy.png
   :name: sod_shock_energy
   :alt:
   :width: 60%


Toro test 2
-----------

Toro's test problem #2 tests robustness to diverging flows and near-zero densities.
The initial state has constant density and temperature, but a jump in velocity.
Left state :math:`\left(\rho_L, u_L, p_L\right) = \left(1.0, -2.0, 0.4\right)` Right
state :math:`\left(\rho_R, u_R, p_R\right) = \left(1.0, 2.0, 0.4\right)`. The result
in a domain of length 5 at time :math:`t=0.6` is shown below.

.. figure:: figs/toro-2.png
   :name: toro-2
   :alt:
   :width: 60%


Toro test 3
-----------

Toro's test problem #3 contains a strong shock close to a contact
discontinuity.  Left initial state :math:`\left(\rho_L, u_L, p_L\right) =
\left(1.0, 0, 1000.0\right)` Right state :math:`\left(\rho_R, u_R,
p_R\right) = \left(1.0, 0, 0.01\right)`.  Time :math:`t = 0.04`.

When evolving pressure, the density peak does not converge to the
analytic solution (solid black line):

.. figure:: figs/toro-3.png
   :name: toro-3
   :alt:
   :width: 60%

However by evolving energy the result converges towards the analytic
solution:

.. figure:: figs/toro-3-energy.png
   :name: toro-3-energy
   :alt:
   :width: 60%

Toro test 4
-----------

Toro's test problem #4 produces two right-going shocks with a contact
between them.  Left state :math:`\left(\rho_L, u_L, p_L\right) =
\left(5.99924, 19.5975, 460.894\right)` Right state
:math:`\left(\rho_R, u_R, p_R\right) = \left(5.99242, -6.19633,
46.0950\right)`.  Result at time :math:`t = 0.15`.

Toro test 5
-----------

The initial conditions for Toro's test problem #5 are the same as test
#3, but the whole system is moving to the left at a uniform speed. The
velocity is chosen so that the contact discontinuity remains almost
stationary at the initial jump location.  Left state
:math:`\left(\rho_L, u_L, p_L\right) = \left(1, 19.59745,
1000.0\right)` Right state :math:`\left(\rho_R, u_R, p_R\right) =
\left(1, -19.59745, 0.01\right)`.  Result at time :math:`t = 0.03`.
