.. _sec-transport_2d

Tokamak axisymmetric transport
==============================

Simulations of transport in axisymmetric tokamak geometries, with
cross-field diffusion and interaction of plasma with neutral gas.

Finding steady state solutions
------------------------------

These models can be run as a time-dependent problem, for example to
study power transients, but the primary application is to finding
steady-state solutions.

Backward Euler solver
~~~~~~~~~~~~~~~~~~~~~

This solver uses PETSc to solve the nonlinear system of equations,
with a Backward Euler timestep to improve the condition number. There
are many choices of algorithm and settings, so the following are
guidelines and may not be optimal for all cases.

.. code-block:: ini

   [solver]
   type = beuler           # Backward Euler steady-state solver
   snes_type = newtonls    # Nonlinear solver
   ksp_type = gmres        # Linear solver
   max_nonlinear_iterations = 10
   pc_type = hypre         # Preconditioner type
   pc_hypre_type = euclid  # Hypre preconditioner type
   lag_jacobian = 500      # Iterations between jacobian recalculations
   atol = 1e-7             # Absolute tolerance
   rtol = 1e-5             # Relative tolerance

PETSc can print quite extensive performance diagnostics. These can be enabled
by putting in the BOUT.inp options file:

.. code-block:: ini

   [petsc]
   log_view = true

This section can also be used to set other PETSc flags, just omitting
the leading `-` from the PETSc option.

   
cvode solver
~~~~~~~~~~~~

CVODE is primarily intended for high-accuracy time integration, rather
than finding steady-state solutions, but can be effective and quite
robust. It tends to struggle at high order, so here we limit it to a
maximum of 3rd order:

.. code-block:: ini

   [solver]
   type = cvode
   use_precon = true   # Use the user-provided preconditioner
   mxstep = 1e5
   mxorder = 3         # Limit to 3rd order
   atol = 1e-12
   rtol = 1e-5

Here `use_precon = true` tells the solver to use the Hermes-3
preconditioners, which are implemented in some components. This
includes preconditioning of parallel heat conduction, and of
cross-field diffusion of neutrals.


Mesh interpolation
~~~~~~~~~~~~~~~~~~

A useful strategy is to start with a low resolution grid, run until
close to steady-state, then interpolate the solution onto a finer mesh
and restart. This process can be repeated as a kind of simplified
multigrid method.



Post-processing
---------------

