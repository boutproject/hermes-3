.. _sec-numerics_performance:

Solvers, performance and numerical methods
==========================================


Solvers
------------------------------

While Hermes-3 has access to any solver within BOUT++ (`see the long
list in the documentation <https://
bout-dev.readthedocs.io/en/stable/user_docs/time_integration.html>`__),
the development and optimisation efforts focuses on two solvers:

CVODE
   A high-order time accurate solver. CVODE capability was initially
   developed for 3D turbulence simulations. By default, it's configured
   to be matrix-free and preconditioned using physical preconditioners
   for electron conduction and neutral cross-field diffusion. It has
   an adaptive order algorithm and will try to go as high order
   in time as possible and take long timesteps. The timestep is limited either
   by the truncation error (the time integration error) or, more commonly, by solver failures.
   See the `developer page <https://computing.llnl.gov/
   projects/sundials/cvode>`__ for details. CVODE is required for 3D
   turbulence simulations. For 2D or 1D, it is generally very slow, but
   tends to be more robust than SNES depending on settings. `The BOUT++ implementation is here. <https
   ://github.com/mikekryjak/BOUT-dev/blob/master/src/solver/impls/cvode/cvode.cxx>`__

SNES
   This an implementation of a first-order in time backward Euler method using the
   `SNES <https://petsc.org/release/manual/snes/>`__
   nonlinear solver in `PETSc <https://petsc.org/release/>`__. It features a range
   of linear solvers and preconditioners available through PETSc, and is typically
   set up with an LU/ILU preconditioner like Hypre ilu, MUMPS or STRUMPACK.
   This allows the entire system to be preconditioned numerically, instead of
   relying on physical preconditioning of select terms as in CVODE.
   Unlike CVODE, the SNES solver features no time integration error control,
   and its discretisation in time is first order only. This is to prioritise
   faster steady-state convergence at the cost of accuracy in the time evolution.

   SNES can lead to substantially improved performance in 1D or 2D
   systems. It is not suitable for 3D turbulence applications due to
   ILU/LU not scaling well to systems that size and due to the vorticity equation
   leading to a very dense Jacobian.

   `The BOUT++ implementation is here. <https://github.com/mikekryjak/
   BOUT-dev/blob/master/src/solver/impls/snes/snes.cxx>`__. See the
   higher level BOUT++ PETSc implementation `here <https://github.com
   /mikekryjak/BOUT-dev/blob/master/src/solver/impls/petsc/petsc.cxx>`__.

Both CVODE and SNES have an adaptive timestepper and will take as long a
timestep as the solver failure rate will allow.

Configuring CVODE
~~~~~~~~~~~~~~~~~~~~~


.. code-block:: ini

   [solver]
   type = cvode
   cvode_precon_method = user  # Use the user-provided preconditioner
   mxstep = 1e9                # Prevent timeout
   mxorder = 3                 # Limit to 3rd order
   atol = 1e-7
   rtol = 1e-5

cvode_precon_method = user
   Enables the physics-based preconditioners. The default is ``none``; other
   choices are ``auto``, ``petsc`` and ``bbd``. The old ``use_precon`` option
   is no longer accepted by CVODE and will raise an error if set.

mxstep = 1e9
   Prevents the solver from timing out if it takes too long to converge.
   1e9 is the maximum number of nonlinear steps allowed for this setting.

mxorder = 3
   Limits the maximum order of the solver to 3rd order. This is useful
   for preventing the solver from getting "stuck" at high order with a
   small timestep. You can alternatively set ``stablimdet = true`` to
   enable automatic `stability limit detection <https://sundials.readthedocs
   .io/en/latest/cvode/Mathematics_link.html#cvode-mathematics-stablimit>`__
   to prevent this.

atol = 1e-7
   Absolute tolerance for the solver. This is the maximum absolute error
   allowed and is required to ensure accuracy/stability when variables are small.
   Available in all solvers.

rtol = 1e-5
   Relative tolerance for the solver, the primary driver of convergence.
   Available in all solvers.


Configuring SNES
~~~~~~~~~~~~~~~~~~~~~

Here is the recommended configuration for the ``SNES`` solver using MUMPS
LU preconditioning (recipe SNES-MUMPS-3):

.. code-block:: ini

   [solver]
   type = snes
   snes_type = newtonls
   ksp_type = fgmres
   pc_type = lu
   max_nonlinear_iterations = 10
   atol = 1e-3
   rtol = 1e-5
   stol = 1e-12
   maxf = 2000
   maxl = 260
   matrix_free_operator = false
   lag_jacobian = 3
   timestep = 0.1
   max_snes_failures = 20
   scale_vars = true

   timestep_control = pid_nonlinear_its
   target_its = 7
   kP = 0.233
   kI = 0.133
   kD = 0.0

   [petsc]
   snes_fd_color_use_mat
   ksp_type = fgmres
   pc_type = lu
   pc_factor_mat_solver_type = mumps
   snes_linesearch_type = basic
   mat_mumps_icntl_24 = 1


Tolerances
   In SNES, ``atol``, ``rtol`` and ``stol`` are the absolute, relative and step nonlinear solver tolerances, respectively.
   This is in contrast to CVODE which has no ``stol`` and where ``atol`` and ``rtol`` are tolerances for time integration accuracy
   with the nonlinear tolerances being derived from them.
   ``stol`` allows a Newton step to succeed if the system has changed by only a very small amount during the step.
   It should be used with caution, as an ``stol`` that is too loose can lead to the solver prematurely declaring convergence,
   leading to an unstable solution.

Jacobian lagging
   One of the largest drivers of performance is the lagging of the Jacobian, which is a linearisation of the nonlinear system
   of equations calculated by finite differences. This is expensive and can cost ~70% of the total computational time within
   a typical timestep.
   When `lag_jacobian = 0`, the Jacobian is updated at every nonlinear iteration, and higher numbers
   indicate the number of nonlinear iterations for which the Jacobian is reused before being recalculated.
   Lower values lead to a higher cost per iteration, but a more accurate linear system, and the optimal
   value is problem dependent. Note that the Jacobian is always recalculated after a nonlinear solver failure.

Adaptive timestep control
   SNES currently uses a PID-based adaptive timestep control strategy implemented on the BOUT++ side.
   By default, it controls the timestep based on the number of nonlinear iterations required
   for convergence as a proxy of convergence difficulty. The three PID coefficients
   ``kP``, ``kI`` and ``kD`` control the proportional, integral and derivative contributions to the timestep adjustment.


PETSc can print quite extensive performance diagnostics. These can be enabled
by putting in the BOUT.inp options file:

.. code-block:: ini

   [petsc]
   log_view = true

This section can also be used to set other PETSc flags, just omitting
the leading ``-`` from the PETSc option.


Performance
------------------------------

Troubleshooting performance
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Hermes-3 uses fully implicit Newton-Krylov solvers. These have high performance
potential but are sensitive to:

- Discontinuities in reaction rates, slope limiter functions, floors or boundary conditions
- Numerical oscillations
- Poor preconditioning

First establish whether the bottleneck is in the linear or the nonlinear solver.
Many linear iterations per nonlinear iteration usually indicate poor preconditioning.
Many nonlinear iterations per timestep, or frequent nonlinear failures, indicate a
strongly nonlinear or stiff system. In CVODE, there are also truncation error failures,
which indicate that the timestep was too large for the desired time evolution accuracy.
Note that this can be triggered both by genuinely fast evolving system, but also by
numerical modes or unphysical fluctuations.

To count iterations:

- Every solver prints a progress line at each output, to the screen and to
  ``BOUT.log.0``, including the number of right-hand-side (RHS) evaluations since
  the previous output.
- With SNES, ``diagnose = true`` in ``[solver]`` also prints a line per timestep
  with the timestep, the nonlinear (``nl iter``) and linear (``lin iter``)
  iterations, the PETSc convergence reason and the number of failures.
- With CVODE, the output file holds cumulative counts: ``cvode_nsteps``
  (timesteps), ``cvode_nniters`` (nonlinear iterations), ``cvode_nliters`` (linear
  iterations), ``cvode_nfevals`` (RHS evaluations), and ``cvode_num_fails`` and
  ``cvode_nonlin_fails`` (truncation error and nonlinear solver failures, respectively).
  Their differences between outputs give the counts per output interval.

PETSc can print the nonlinear residual norm at each iteration. This is the scalar quantity
that the solver tries to drive to zero at each iteration:

.. code-block:: ini

   [petsc]
   snes_monitor
   snes_converged_reason

``snes_monitor`` prints the 2-norm of the residual at every nonlinear iteration, and
``snes_converged_reason`` prints why each solve stopped. Close to convergence, the
norm should fall faster with each iteration. A norm that stalls or rises suggests an
inaccurate Jacobian (for example from ``lag_jacobian`` set too high) or a
non-differentiable term such as a floor, limiter or switch. ``log_view`` (see
Configuring SNES) shows where the time is spent.

The nonlinear residual also shows which equations and grid cells cause the failures.
Set ``diagnose = true`` in the ``[solver]`` section to write the residual of each
evolved variable at every output time as ``resid_<variable>`` (e.g. ``resid_Pe`` or
``resid_Nd+``). This is the residual that SNES drives below ``atol`` and ``rtol``,
taken from the last solve before the output. Large values mark the equations and
cells that dominate the error. Setting ``diagnose_failures = true`` in ``[solver]``
also prints the range of each variable and of its time derivative after every
failed solve.

Neutral model performance
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The neutral model can cause significant performance degradation when using CVODE due to the poor performance
of the current physical preconditioner.
This is likely due to the short transport timescales in the perpendicular transport terms when crossing the separatrix,
with faster fluxes generally causing more trouble for the solver.

One way to mitigate this is to reduce the ``flux_limit`` value to artificially suppress neutral diffusive fluxes.
However, this can significantly affect the physical accuracy of the solution.

Significant performance improvements have been achieved in both CVODE and SNES by lagging a part of the
neutral flux limiter formula:

.. code-block:: ini

   [d]
   type = neutral_mixed
   lag_limiter_gradient = true

This computes the regularised gradients in the limiter denominators (see
:ref:`sec-neutral_mixed`) once per timestep and holds them fixed for linear and
nonlinear iterations within the timestep. In 2D tests, this had a very small impact on
the time evolution of the solution but led to a large performance improvement.
The benefit is problem dependent. The option is off by default.


Numerics
------------------------------

Slope (flux) limiters
~~~~~~~~~~~~~~~~~~~~~

The choice of slope limiter is important: a dissipative limiter increases
numerical dissipation but can substantially improve solver stability and
robustness. See :ref:`sec-slope-limiter-settings` for how to change the
limiter. The default is ``MC``, which has a good balance between
stability and accuracy.

Dynamics parallel to the magnetic field are solved using a 2nd-order
slope-limiter method.  For any number of fluids we solve the number
density :math:`n`, momentum along the magnetic field,
:math:`mnv_{||}`, and either pressure :math:`p` or energy
:math:`\mathcal{E}`. Here :math:`m` is the particle mass, so that :math:`mn`
is the mass density. :math:`v_{||}` is the component of the flow
velocity in the direction of the magnetic field, and is aligned with
one of the mesh coordinate directions.  All quantities are cell
centered.

Cell edge values are reconstructed using a configurable limiter. The
available compile-time choices of ``HERMES_SLOPE_LIMITER`` are:

``MC``
   Monotonised central second-order limiter. This is the default.

``VanAlbada``
   Symmetric second-order limiter using a smooth approximation to the
   usual limiter switch. This is often friendlier to nonlinear solvers
   and finite-difference Jacobian calculations because the reconstruction
   remains differentiable at extrema.

``WENO3``
   Third-order Jiang-Shu WENO reconstruction on the same three-point stencil.
   This is generally smoother and less dissipative than TVD limiters, but it
   does not enforce strict monotonicity.

``MinMod``
   More dissipative second-order TVD limiter.

``Upwind``
   First-order upwind reconstruction. This is the most dissipative, but is
   sometimes useful for robustness tests.

``Superbee``
   Aggressive second-order TVD limiter with reduced diffusion in sharp fronts,
   but it can be less robust.

If :math:`f_i` is the value of field :math:`f` at the
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

.. _sec-slope-limiter-settings:

Selecting slope limiter and conduction method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Hermes-3 selects the parallel advection slope limiter at compile time
through CMake cache options:

.. code-block:: bash

   cmake -S . -B build \
     -DHERMES_SLOPE_LIMITER=MC

``HERMES_SLOPE_LIMITER``
   Available values are ``MC``, ``VanAlbada``, ``WENO3``, ``MinMod``,
   ``Upwind``, and ``Superbee``.


The
``Div_par_K_Grad_par_mod`` conduction discretisation can be choosen at run time with:
``conduction_method``
   Selects the parallel heat conduction discretisation used in
   :cpp:func:`Div_par_K_Grad_par_mod`. Available values are:

   ``Original``
      Separately averages :math:`K`, :math:`J`, :math:`g_{22}`, and
      :math:`dy` at the face and then multiplies them together.

   ``ProductJK``
      Uses the same stencil as ``Original`` but averages :math:`J K`
      together before applying the face gradient.

   ``Harmonic``
      Uses a harmonic average of the half-cell conductances
      :math:`K J / (g_{22} dy)`. This better matches a series-resistance
      interpretation of two adjacent half-cells and can give noticeably
      different results when coefficients or cell sizes vary strongly.

The selected values are reported at startup in the Hermes-3 log and are
also recorded in the options tree as ``slope_limiter`` and
``conduction_method``.

Controlling Lax flux strength with sound_speed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _sec-sound_speed:

**sound_speed**

By default, the Lax flux strength is calculated based on the sound speed
of each species individually. Instead, the ``sound_speed`` component can be used
component to sum the pressure of all species and divide by the sum of the mass density
of all species:

.. math::

   c_s = \sqrt{\sum_i P_i / \sum_i m_in_i}

This is set in the state as `sound_speed`. It is highly recommended to use
this when evolving electron momentum, as the electrons will represent the fastest
sound speed in the system, increasing the amount of numerical damping from the Lax flux.
This is enabled by default.

NOTE:
   When using this component **without** electron momentum being evolved,
   set ``electron_dynamics = false`` in the ``[sound_speed]`` section to prevent
   the electron timescale from being applied, or you will have too much damping.

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

.. _sec-low-n-sources:

Low-density sources (experimental)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``low_n_sources`` component is an **experimental feature** intended for use
in simulations where low species densities are expected to cause numerical
issues. It adds artificial sources and sinks via ``pseudo-reactions`` to restore
depleted species towards a density threshold. Note that, while these terms are
implemented via the regular :ref:`reactions framework<sec-reactions-dev>`, they
are distinct from anything in the :code:`[reactions]` block, and are included
purely to improve numerical stability.

To enable the sources, add :code:`low_n_sources` to the component list in your
input file:

.. code-block:: ini

   [hermes]
   components = ..., low_n_sources

   [low_n_sources]
   floor_factor = 10
   diagnose = true

The activation threshold for each source is the :code:`density_floor` of its target
species multiplied by a :code:`floor_factor` (default :code:`10`). Below this density,
the pseudo-reaction adds density sources and sinks to the participating species,
together with the associated momentum and energy transfers. They do not impose a
hard lower bound on density. Setting :code:`diagnose = true` adds the density
sources associated with each pseudo-reaction to the simulation output.

Pseudo-reactions
^^^^^^^^^^^^^^^^

The types of pseudo-reactions available in each simulation are determined
automatically based on which species are present. Reactants and products are
typically, but not necessarily, the same as the equivalent physical reactions.
The pseudo-reaction types are listed below, with example reaction strings in
parentheses for a simulation including :math:`H`, :math:`H^+`, :math:`H_2` and :math:`H_2^+`.

* **Recombination** (:math:`H^+ + e \to H` and :math:`H_2^+ + e \to H_2`):

  * Restores depleted neutral species, consuming some ions.
  * ``Available when`` at least one neutral and singly charged ion pair is present.
  * ``Triggered when`` the neutral (atom or molecule) population falls below its threshold.
  * ``Boosted when`` the product species is an atom and the corresponding molecular species (if present) is also depleted.

* **Ionisation** (:math:`H \to H^+ + e` and :math:`H_2 \to H_2^+ + e`):

  * Restores depleted ion species, consuming some neutrals.
  * ``Available when`` at least one pair of neutral and ion species are present in the simulation.
  * ``Triggered when`` the ion population falls below its threshold.
  * Smoothly ``suppressed when`` the reactant species is below its density threshold.

* **Association** (:math:`2H \to H_2`):

  * Restores depleted molecules, consuming some atoms.
  * ``Available when`` at least one pair of atom and molecule species are present in the simulation.
  * ``Triggered when`` the molecular population falls below its threshold.

* **External ion source** (represented internally as :math:`H^+ \to 2H^+` and :math:`H_2^+ \to 2H_2^+`):

  * Adds ions to the system (introduces a net source of nuclei and mass).
  * ``Available when`` at least one ion species is present in the simulation.
  * ``Triggered only as a last resort when`` the ion and corresponding neutral populations are simultaneously below their thresholds.
  * If the corresponding neutral species is not present, the external source is triggered solely by ion deficit.
  * Adds momentum and thermal energy at the local ion velocity and temperature. These terms can be disabled with :code:`external_src_mom_energy_enabled = false`.
  * ``Boosted when`` a corresponding neutral species is included in the simulation and is also depleted.

Note that, for now, only singly charged ions can be included in pseudo-reactions, other charge states are ignored.

Choosing the restoration timescale
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :code:`timescale` option controls how rapidly pseudo-reactions restore
species densities towards their threshold(s). User-supplied values must be
in normalised code units; the default corresponds to :math:`10^{-5}` s.
When setting this option, aim for a timescale short enough that
pseudo-reactions are able to compete with physical reactions in the low-density
regime where they are triggered, but not so short that they behave like a hard
density floor, which may increase stiffness and worsen solver convergence.
