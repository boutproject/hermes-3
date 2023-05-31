.. _sec-getting_started:

Getting started
===============

Installing
----------

Only CMake is supported for building Hermes-3 and running the tests.
During configuration `BOUT++
<https://github.com/boutproject/BOUT-dev/>`_ will be automatically
downloaded as a submodule, together with some dependencies. `NetCDF
<https://www.unidata.ucar.edu/software/netcdf/>`_ and `FFTW
<https://www.fftw.org/>`_ are assumed to be installed already.  The
`SUNDIALS <https://computing.llnl.gov/projects/sundials>`_ library is
strongly recommended for time-dependent simulations, and `PETSc
<https://petsc.org>`_ is needed to run some of the steady-state
transport solver examples.

If you only want to run time-dependent simulations, then the
recommended way to build Hermes-3 links to the `SUNDIALS
<https://computing.llnl.gov/projects/sundials>`_ library:

#. Configure with cmake, downloading and linking to SUNDIALS:

   .. code-block:: bash

      cmake . -B build -DBOUT_DOWNLOAD_SUNDIALS=ON

#. Build, compiling Hermes-3 and all dependencies:

   .. code-block:: bash

      cmake --build build

#. Run the unit and integrated tests to check that everything is working:

   .. code-block:: bash

      cd build
      ctest

Note that the integrated tests require MPI, and so may not run on the
head nodes of many computing clusters.

The CMake configuration can be customised: See the [BOUT++
documentation](https://bout-dev.readthedocs.io/en/latest/user_docs/installing.html#cmake)
for examples of using `cmake` arguments, or edit the compile options
interactively before building:

.. code-block:: bash

   ccmake . -B build

If you have already installed BOUT++ and want to use that rather than
configure and build BOUT++ again, set `HERMES_BUILD_BOUT` to `OFF` and pass
CMake the path to the BOUT++ `build` directory e.g.

.. code-block:: bash

   cmake . -B build -DHERMES_BUILD_BOUT=OFF -DCMAKE_PREFIX_PATH=$HOME/BOUT-dev/build

Note that Hermes-3 currently requires BOUT++ version 5.

Building with PETSC
-------------------

When building PETSc it is recommended to include ``hypre``. The
following PETSc configure line is a good starting point:

.. code-block:: bash

   ./configure --with-mpi=yes --download-hypre --download-make --with-fortran-bindings=0

To configure Hermes-3 with PETSc, use the ``-DBOUT_USE_PETSC=ON`` flag:

.. code-block:: bash

      cmake . -B build -DBOUT_DOWNLOAD_SUNDIALS=ON -DBOUT_USE_PETSC=ON

If the ``PETSC_DIR`` and ``PETSC_ARCH`` environment variables have been set,
then CMake should pick them up.

Numerical methods
-----------------

Advection operators in Hermes-3 use slope limiters, also called `flux
limiters <https://en.wikipedia.org/wiki/Flux_limiter` to suppress
spurious numerical oscillations near sharp features, while converging
at 2nd-order in smooth regions. In general there is a trade-off
between suppression of numerical oscillations and dissipation: Too
little dissipation results in oscillations that can cause problems
(e.g. negative densities), while too much dissipation smooths out real
features and requires higher resolution to converge to the same
accuracy. The optimal choice of method is problem-dependent.

The CMake option ``HERMES_SLOPE_LIMITER`` sets the choice of slope
limiter.  The default method is ``MinMod``, which has been found to
provide a good balance for problems of interest. If less dissipation
is required then this can be changed to ``MC`` (for Monotonized
Central); For more dissipation (but 1st-order convergence) change it
to ``Upwind``.
