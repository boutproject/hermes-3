.. _sec-numerical-methods:

Numerical methods
=================

Parallel dynamics
-----------------

Dynamics parallel to the magnetic field are solved using a 2nd-order
slope-limiter method.  For any number of fluids we solve the number
density :math:`n`, momentum along the magnetic field,
:math:`mnv_{||}`, and either pressure :math:`p` or energy
:math:`\mathcal{E}`. Here :math:`m` is the particle mass, so that $mn$
is the mass density. :math:`v_{||}` is the component of the flow
velocity in the direction of the magnetic field, and is aligned with
one of the mesh coordinate directions.  All quantities are cell
centered.

Cell edge values are by default reconstructed using a MinMod method
(other limiters are available, including 1st-order upwind, Monotonized
Central, and Superbee). If :math:`f_i` is the value of field :math:`f` at the
center of cell :math:`i`, then using MinMod slope limiter the gradient :math:`g_i`
inside the cell is:

.. math::

   g_i = \left\{\begin{array}{ll}
   0 & \textrm{if $\left(f_{i+1} - f_{i}\right) \left(f_{i} - f_{i-1}\right) < 0$} \\
   f_{i+1} - f_{i} & \textrm{if $\left|f_{i+1} - f_{i}\right| < \left|f_{i} - f_{i-1}\right|$} \\
   f_{i} - f_{i-1} & \textrm{Otherwise}
   \end{array}\right.

The values at the left and right of cell :math:`i` are:

.. math::

   \begin{align}
   f_{i, R} &= f_i + g_i / 2 \nonumber \\
   f_{i, L} &= f_i - g_i / 2
   \end{align}

This same reconstruction is performed for :math:`n`, :math:`v_{||}` and :math:`p` (or
:math:`\mathcal{E}`). The flux :math:`\Gamma_{i+1/2}` between cell :math:`i` and :math:`i+1`
is:

.. math::

   \Gamma_{f, i+1/2} = \frac{1}{2}\left(f_{i,R} v_{||i,R} + f_{i+1,L}v_{||i+1,L}\right) + \frac{a_{max,i+1/2}}{2}\left(f_{i,R} - f_{i+1,L}\right)

This includes a Lax flux term that penalises jumps across cell edges,
and depends on the maximum local wave speed, :math:`a_{max}`. Momentum is
not reconstructed at cell edges; Instead the momentum flux is
calculated from the cell edge densities and velocities:

.. math::

   \Gamma_{nv, i+1/2} = \frac{1}{2}\left(n_{i,R} v_{||i,R}^2 + n_{i+1,L}v_{||i+1,L}^2\right) + \frac{a_{max,i+1/2}}{2}\left(n_{i,R}v_{||i,R} - n_{i+1,L}v_{||i+1,R}\right)

The wave speeds, and so :math:`a_{max}`, depend on the model being solved,
so can be customised to e.g include or exclude Alfven waves or
electron thermal speed. For simple neutral fluid simulations it is:

.. math::

   a_{max, i+1/2} = \max\left(\left|v_{||i}\right|, \left|v_{||i+1}\right|, \sqrt{\frac{\gamma p_{i}}{mn_i}}, \sqrt{\frac{\gamma p_{i+1}}{mn_{i+1}}}\right)

The divergence of the flux, and so the rate of change of :math:`f` in cell
:math:`i`, depends on the cell area perpendicular to the flow, :math:`A_i`, and cell volume :math:`V_i`:

.. math::

   \nabla\cdot\left(\mathbf{b} f v_{||}\right)_{i} = \frac{1}{V_i}\left[\frac{A_{i} + A_{i+1}}{2}\Gamma_{f, i+1/2} - \frac{A_{i-1} + A_{i}}{2}\Gamma_{f, i-1/2}\right]

Boundaries
~~~~~~~~~~

At boundaries along the magnetic field the flow of particles and
energy are set by e.g.  Bohm sheath boundary conditions or no-flow
conditions. To ensure that the flux of particles is consistent with
the boundary condition imposed at cell boundaries, fluxes of density
:math:`n` and also :math:`p` or :math:`\mathcal{E}` are set to the simple mid-point
flux:

.. math::

   \Gamma_{f, i+1/2}^{boundary} = f_{i+1/2}v_{||i+1/2}

where :math:`f_{i+1/2} = \frac{1}{2}\left(f_{i} + f_{i+1}\right)` and
:math:`v_{||i+1/2} = \frac{1}{2}\left(v_{||i} + v_{||i+1}\right)` are the
mid-point averages where boundary conditions are imposed.  It has been
found necessary to include dissipation in the momentum flux at the
boundary, to suppress numerical overshoots due to the narrow boundary
layers that can form:

.. math::

   \Gamma_{nv, i+1/2}^{boundary} = n_{i,R}v_{||i,R}v_{||i+1/2} + a_{max}\left[n_{i,R}v_{||i,R} - n_{i+1/2}v_{||i+1/2}\right]

where :math:`n_{i+1/2} = \frac{1}{2}\left(n_{i} + n_{i+1}\right)`.

