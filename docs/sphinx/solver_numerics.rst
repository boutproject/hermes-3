.. _sec-solver_numerics:

Solvers and numerical methods
==============================


Solvers
------------------------------

While Hermes-3 has access to any solver within BOUT++ (`see the long
list in the documentation <https://
bout-dev.readthedocs.io/en/stable/user_docs/time_integration.html>`_), 
the development and optimisation efforts focuses on two solvers:

CVODE
   A high-order time accurate solver. CVODE capability was initially
   developed for 3D turbulence simulations. By default, it's configured
   to be matrix-free and preconditioned using physical preconditioners
   for electron conduction and neutral cross-field diffusion. It has
   an adaptive order algorithm and will try to go as high order
   in time as possible and take long timesteps. Please
   see the `developer page <https://computing.llnl.gov/
   projects/sundials/cvode>`_ for details. CVODE is required for 3D
   turbulence simulations and currently gives better results in 2D 
   simulations. `The BOUT++ implementation is here. <https
   ://github.com/mikekryjak/BOUT-dev/blob/master/src/solver/impls/cvode/cvode.cxx>``

beuler
   A first-order in time implementation of the backward Euler method.
   This is implemented using the `SNES <https://petsc.org/release/manual/snes/>`_ 
   nonlinear solver in `PETSc <https://petsc.org/release/>`_, and 
   currently uses a GMRES linear solver and an ILU-type preconditioner
   in `hypre <https://hypre.readthedocs.io/en/latest/ch-intro.html>`_. 
   The algebraic preconditioner allows for extreme speedups over CVODE,
   and can already be used in 1D with good results. 2D performance is
   highly mixed and still in development. 3D use is not feasible as the 
   matrix size precludes the use of algebraic preconditioners.
   Due to its first order nature, it will be less accurate for time-dependent
   problems, and it is mostly intended to solve for steady-state problems.
   `The BOUT++ implementation is here. <https://github.com/mikekryjak/
   BOUT-dev/blob/master/src/solver/impls/snes/snes.cxx>`_. See the 
   higher level BOUT++ PETSc implementation `here <https://github.com
   /mikekryjak/BOUT-dev/blob/master/src/solver/impls/petsc/petsc.cxx>`_.

Both CVODE and beuler have an adaptive timestepper and will take as long a 
timestep as the Newton iteration failure rate will allow.

Configuring CVODE
~~~~~~~~~~~~~~~~~~~~~


.. code-block:: ini

   [solver]
   type = cvode
   use_precon = true   # Use the user-provided preconditioner
   mxstep = 1e9        # Prevent timeout
   mxorder = 3         # Limit to 3rd order
   atol = 1e-7
   rtol = 1e-5

use_precon = true
   Enables the physics-based preconditioners

mxstep = 1e9
   Prevents the solver from timing out if it takes too long to converge.
   1e9 is the maximum number of nonlinear steps allowed for this setting.

mxorder = 3
   Limits the maximum order of the solver to 3rd order. This is useful
   for preventing the solver from getting "stuck" at high order with a
   small timestep. You can alternatively set ``stablimdet = true`` to
   enable automatic `stability limit detection <https://sundials.readthedocs
   .io/en/latest/cvode/Mathematics_link.html#cvode-mathematics-stablimit>`_ 
   to prevent this.

atol = 1e-7
   Absolute tolerance for the solver. This is the maximum absolute error 
   allowed and is required to ensure accuracy/stability when variables are small.
   Available in all solvers.

rtol = 1e-5
   Relative tolerance for the solver, the primary driver of convergence.
   Available in all solvers.


Configuring beuler
~~~~~~~~~~~~~~~~~~~~~

Here is the re
.. code-block:: ini

   [solver]
   type = beuler                 # Backward Euler steady-state solver
   snes_type = newtonls          # Nonlinear solver
   ksp_type = gmres              # Linear solver
   max_nonlinear_iterations = 10 # Max Newton iterations until iteration failure
   pc_type = hypre               # Preconditioner type
   pc_hypre_type = euclid        # Hypre preconditioner type
   lag_jacobian = 500            # Iterations between jacobian recalculations
   atol = 1e-7                   # Absolute tolerance
   rtol = 1e-5                   # Relative tolerance


There is an additional option ``stol`` which allows the Newton iteration
to "converge" if the simulation has only changed by a small amount. 
This enables the simulation to iterate almost instantly if it detects that the 
results aren't changing. 

PETSc can print quite extensive performance diagnostics. These can be enabled
by putting in the BOUT.inp options file:

.. code-block:: ini

   [petsc]
   log_view = true

This section can also be used to set other PETSc flags, just omitting
the leading `-` from the PETSc option.

   

Numerics
------------------------------

Slope (flux) limiters 
~~~~~~~~~~~~~~~~~~~~~

Dynamics parallel to the magnetic field are solved using a 2nd-order
slope-limiter method.  For any number of fluids we solve the number
density :math:`n`, momentum along the magnetic field,
:math:`mnv_{||}`, and either pressure :math:`p` or energy
:math:`\mathcal{E}`. Here :math:`m` is the particle mass, so that :math:`mn`
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

Controlling Lax flux strength
~~~~~~~~~~

See the ``sound_speed`` component.

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


