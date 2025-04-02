.. _sec-installation_using_cmake:

Installation using CMake
===============

Compilation process
----------

Compilation is achieved in two stages - the first is configuration where all the compile-time
options are read in. The second is the build which results in a ready-to-use Hermes-3 installation
in a directory named `build` by default. The build directory name can be changed to have
multiple builds available at the same time.

If you make changes to the code, you can skip straight to the build stage to save time.
Only modified files will be recompiled.

Hermes-3 is built using `CMake <https://cmake.org>`_. During configuration `BOUT++
<https://github.com/boutproject/BOUT-dev/>`_ will be automatically
downloaded as a submodule, together with some dependencies. The correct version 
of `netCDF <https://www.unidata.ucar.edu/software/netcdf/>`_ is downloaded 
and compiled automatically for convenience. `FFTW
<https://www.fftw.org/>`_ is assumed to be installed already. 

Hermes-3 uses two solvers: `SUNDIALS <https://computing.llnl.gov/projects/sundials>`_ `cvode` for
time-dependent simulations and the faster `PETSc
<https://petsc.org>`_ `beuler` for steady-state transport problems. While SUNDIALS
can be downloaded and installed automatically, PETSc requires manual installation.

Installing with SUNDIALS (cvode) only
----------

If you only want to use the `cvode` solver, then the
recommended way to build Hermes-3 links to the SUNDIALS library:


1. Configure with cmake, downloading and linking to SUNDIALS and NetCDF4:

   .. code-block:: bash

      cmake . -B build -DBOUT_DOWNLOAD_SUNDIALS=ON -DBOUT_DOWNLOAD_NETCDF_CXX4=ON

2. Build, compiling Hermes-3 and all dependencies using 4 parallel cores
(adjust as necessary):

   .. code-block:: bash

      cmake --build build -j 4

3. Run the unit and integrated tests to check that everything is working:

   .. code-block:: bash

      cd build
      ctest

Note that the integrated tests require MPI, and so may not run on the
head nodes of many computing clusters.


Installing with SUNDIALS and PETSc (beuler)
-------------------

The steady-state solver beuler requires PETSc and is often preconditioned using the `hypre`
package, which is automatically downloaded and configured during PETSc installation.

Here is an example PETSc configure setup:

.. code-block:: bash

   ./configure --with-mpi=yes --download-hypre --download-make --with-fortran-bindings=0 --with-debugging=0

Here is an example working script to automatically download and compile PETSc on `Viking2`:

.. code-block:: bash

      mkdir petsc-build
      wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.17.4.tar.gz
      tar xzf petsc-3.17.4.tar.gz
      cd petsc-3.17.4
      ./configure COPTFLAGS="-O3" CXXOPTFLAGS="-O3" FOPTFLAGS="-O3" --download-hypre --with-debugging=0 --prefix=../petsc-build
      make -j 4 PETSC_DIR=$PWD PETSC_ARCH=arch-linux-c-opt all
      make -j 4 PETSC_DIR=$PWD PETSC_ARCH=arch-linux-c-opt install
      make -j 4 PETSC_DIR=$PWD/../petsc-build PETSC_ARCH="" check

and on `ARCHER2`:

.. code-block:: bash

      mkdir petsc-build
      wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.17.4.tar.gz
      tar xzf petsc-3.17.4.tar.gz
      cd petsc-3.17.4
      ./configure --CC=cc --CXX=CC --FC=ftn COPTFLAGS="-Ofast" CXXOPTFLAGS="-Ofast" FOPTFLAGS="-Ofast" --with-batch --known-64-bit-blas-indices=0 --known-sdor-returns-double=0 --known-snrm2-returns-double=0 --with-fortran-bindings=0 --download-hypre --with-debugging=0 --prefix=../petsc-build
      make -j 4 PETSC_DIR=$PWD PETSC_ARCH=arch-linux-c-opt all
      make -j 4 PETSC_DIR=$PWD PETSC_ARCH=arch-linux-c-opt install
      make -j 4 PETSC_DIR=$PWD/../petsc-build PETSC_ARCH="" check

And here is a working configure example for `Perlmutter`:

.. code-block:: bash

    ./configure \
      --with-mpi=yes --with-precision=double --with-scalar-type=real --with-shared-libraries=1 \
      --with-debugging=0 {C,CXX,F}OPTFLAGS="-O3 -march=native" \
      --download-hypre --download-fblaslapack=1 \
      --prefix=$HOME/local/petsc-3.22.3

Once PETSc is installed, link it to Hermes-3 using the ``-DBOUT_USE_PETSC=ON`` CMake flag:

.. code-block:: bash

      cmake . -B build -DBOUT_DOWNLOAD_SUNDIALS=ON -DBOUT_DOWNLOAD_NETCDF_CXX4=ON -DBOUT_USE_PETSC=ON

If the ``PETSC_DIR`` and ``PETSC_ARCH`` environment variables have been set,
then CMake should pick them up. If it doesn't, try doing a clean build by removing
any previously generated build directories.


Dependencies
-------------------
Since Hermes-3 heavily relies on BOUT++, the `BOUT++ documentation on installation and
dependencies <https://bout-dev.readthedocs.io/en/stable/user_docs/quickstart.html#prerequisites>`_ 
contains a lot of useful information. Below is a selection of working module lists
for several HPC systems. It is recommended you start with a clean module environment 
by executing `module purge` first.

YPI Workstations:

.. code-block:: bash

   module load mpi/OpenMPI/4.1.1-GCC-10.3.0
   module load devel/CMake/3.20.1-GCCcore-10.3.0
   module load numlib/OpenBLAS/0.3.15-GCC-10.3.0
   module load lib/FlexiBLAS/3.0.4-GCC-10.3.0

ARCHER2:

.. code-block:: bash

   module swap PrgEnv-cray/8.3.3
   module swap cce/15.0.0
   module swap cray-mpich/8.1.23
   module load cray-python/3.9.13.1 
   module load netcdf4 
   module load cmake 
   module load cray-hdf5 
   module load cray-netcdf/4.9.0.1 
   module load cray-parallel-netcdf/1.12.3.1 
   module load cray-fftw/3.3.10.3 
   module load valgrind4hpc

Marconi:

.. code-block:: bash

   module load tools/git/2.32.0-GCCcore-10.3.0-nodocs
   module load mpi/OpenMPI/4.1.1-GCC-10.3.0
   module load devel/CMake/3.20.1-GCCcore-10.3.0
   module load numlib/OpenBLAS/0.3.15-GCC-10.3.0
   module load data/netCDF/4.8.0-gompi-2021a
   module load lang/SciPy-bundle/2021.05-foss-2021a

Viking2:

.. code-block:: bash

   module load OpenMPI/4.1.1-GCC-10.3.0
   module load git/2.32.0-GCCcore-10.3.0-nodocs
   module load CMake/3.20.1-GCCcore-10.3.0
   module load OpenBLAS/0.3.15-GCC-10.3.0
   module load netCDF/4.8.0-gompi-2021a
   module load SciPy-bundle/2021.05-foss-2021a

Ancalagon:

.. code-block:: bash

   module load OpenMPI/4.1.1-GCC-10.3.0 
   module load CMake/3.20.1-GCCcore-10.3.0 
   module load OpenBLAS/0.3.15-GCC-10.3.0 
   module load SciPy-bundle/2021.05-foss-2021a 
   module load netCDF/4.8.0-gompi-2021a

Perlmutter:

.. code-block:: bash

   source /opt/cray/pe/cpe/23.03/restore_lmod_system_defaults.sh
   module load craype-x86-rome
   module load libfabric
   module load craype-network-ofi
   module load xpmem
   module load cray-libsci
   module load PrgEnv-gnu
   module load cray-mpich
   module load python
   module load cray-fftw
   module load cray-hdf5
   module load cray-netcdf


Slope (flux) limiter settings
-----------------

Advection operators in Hermes-3 use slope limiters, also called `flux
limiters <https://en.wikipedia.org/wiki/Flux_limiter>`_ to suppress
spurious numerical oscillations near sharp features, while converging
at 2nd-order in smooth regions. In general there is a trade-off
between suppression of numerical oscillations and dissipation: Too
little dissipation results in oscillations that can cause problems
(e.g. negative densities), while too much dissipation smooths out real
features and requires higher resolution to converge to the same
accuracy. The optimal choice of method is problem-dependent.

The CMake flag ``-DHERMES_SLOPE_LIMITER`` sets the choice of slope
limiter.  The default method is ``MC``, which has been found to
provide a good balance for problems of interest. If more dissipation
is required then this can be changed to ``MinMod``; 
if less dissipation is required then this can be changed
to ``Superbee``.

The appropriate limiter is problem-dependent. ``MinMod`` can work well
for 1D tokamak simulations with steep gradients, e.g. simulations of detachment
transients in high power machines which are already under-dissipative
due to the lack of cross-field transport. The use of ``MinMod`` in 2D or 3D can
lead to over-dissipation, but greater robustness.


Compiling in debug mode
-----------------
Please see the `relevant page <https://bout-dev.readthedocs.io/en/stable/user_docs/advanced_install.html#optimisation-and-run-time-checking>`_ 
in the BOUT++ documentation.


Custom versions of BOUT++
-----------------

If you have already installed BOUT++ and want to use that rather than
configure and build BOUT++ again, set ```-HERMES_BUILD_BOUT=OFF``` and pass
CMake the path to the BOUT++ `build` directory e.g.

.. code-block:: bash

   cmake . -B build -DHERMES_BUILD_BOUT=OFF -DCMAKE_PREFIX_PATH=$HOME/BOUT-dev/build

The version of BOUT++ required by Hermes-3 is periodically updated, and is usually derived 
from a commit on the `next` branch of BOUT++. The up to date commit can be found in the 
`"external" directory of the Hermes-3 repo 
<https://github.com/bendudson/hermes-3/tree/master/external>`_.


Custom configuration of CMake
-----------------

The CMake configuration can be customised: See the `BOUT++
documentation
<https://bout-dev.readthedocs.io/en/latest/user_docs/installing.html#cmake>`_
for examples of using `cmake` arguments, or edit the compile options
interactively before building:

.. code-block:: bash

   ccmake . -B build


Troubleshooting issues
-----------------

The first step to troubleshooting compilation issues should always to delete
build folder for a fresh compilation. This can resolve several types of issues.

There have also been several reported issues due to Conda (e.g. making 
BOUT++ pick up the Conda MPI installation instead of the module one). A 
workaround is to compile with the CMake flag `-DBOUT_IGNORE_CONDA_ENV=ON`.



