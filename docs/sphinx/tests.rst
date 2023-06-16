.. _sec-tests:

Tests
=====

The specification of the Toro tests used here is taken from
`Walker (2012) <https://doi.org/10.1371/journal.pone.0039999>`_,
originally from Toro's book `Riemann Solvers and Numerical Methods for
Fluid Dynamics <https://link.springer.com/book/10.1007/b79761>`_.

1D fluid (MMS)
--------------

``tests/integrated/1D-fluid``

This convergence test using the Method of Manufactured Solutions (MMS)
solves fluid equations in the pressure form:

.. math::

   \begin{aligned}
   \frac{\partial n}{\partial t} &= -\nabla\cdot\left(n\mathbf{b}v_{||}\right) \\
   \frac{\partial p}{\partial t} &= -\nabla\cdot\left(p\mathbf{b}v_{||}\right) - \frac{2}{3}p\nabla\cdot\left(\mathbf{b}v_{||}\right) \\
   \frac{\partial}{\partial t}\left(mnv_{||}\right) &= -\nabla\cdot\left(nv_{||}\mathbf{b}v_{||}\right) - \partial_{||}p
   \end{aligned}


.. figure:: figs/fluid_norm.png
   :name: fluid_norm
   :alt:
   :width: 60%

Sod shock
---------

``tests/integrated/sod-shock`` and ``tests/integrated/sod-shock-energy``

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

Toro test 1
-----------

``tests/integrated/toro-1``

Toro's test problem #1, from `Riemann Solvers and Numerical Methods
for Fluid Dynamics <https://link.springer.com/book/10.1007/b79761>`_
is a variation of Sod's shock tube problem. The left state is moving
into the right, increasing the speed of the resulting shock. Left
state :math:`\left(\rho_L, u_L, p_L\right) = \left(1.0, 0.75,
1.0\right)` Right state :math:`\left(\rho_R, u_R, p_R\right) =
\left(0.125, 0, 0.1\right)`. The size of the domain is 5, and
the reference result is given at time :math:`t = 0.8`.

Toro test 2
-----------

``tests/integrated/toro-2`` and ``tests/integrated/toro-2-energy``

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

``tests/integrated/toro-3`` and ``tests/integrated/toro-3-energy``

Toro's test problem #3 contains a strong shock close to a contact
discontinuity.  Left initial state :math:`\left(\rho_L, u_L, p_L\right) =
\left(1.0, 0, 1000.0\right)` Right state :math:`\left(\rho_R, u_R,
p_R\right) = \left(1.0, 0, 0.01\right)`.  Time :math:`t = 0.04`.

When evolving pressure, the simulation is robust but the density peak
does not converge to the analytic solution (solid black line):

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

``tests/integrated/toro-4`` and ``tests/integrated/toro-4-energy``

Toro's test problem #4 produces two right-going shocks with a contact
between them.  Left state :math:`\left(\rho_L, u_L, p_L\right) =
\left(5.99924, 19.5975, 460.894\right)` Right state
:math:`\left(\rho_R, u_R, p_R\right) = \left(5.99242, -6.19633,
46.0950\right)`.  Result at time :math:`t = 0.15`.

Toro test 5
-----------

``tests/integrated/toro-5`` and ``tests/integrated/toro-5-energy``

The initial conditions for Toro's test problem #5 are the same as test
#3, but the whole system is moving to the left at a uniform speed. The
velocity is chosen so that the contact discontinuity remains almost
stationary at the initial jump location.  Left state
:math:`\left(\rho_L, u_L, p_L\right) = \left(1, -19.59745,
1000.0\right)` Right state :math:`\left(\rho_R, u_R, p_R\right) =
\left(1, -19.59745, 0.01\right)`.  Result at time :math:`t = 0.03`.
