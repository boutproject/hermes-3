.. _sec-installation:

Installation
===============

Hermes-3 can be installed using CMake or Spack. Using CMake is a more manual process 
which requires you to provide all of the dependencies yourself, but it has been used
extensively and the documentation provides a module list for several HPC systems.

Spack capability has recently been added and allows an easy way to install the 
dependencies, making it much easier to install Hermes-3 without a module environment,
e.g. on a workstation or laptop.


Using CMake
----------

Compilation process
~~~~~~~~~~

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
~~~~~~~~~~

If you only want to use the `cvode` solver, then the
recommended way to build Hermes-3 links to the SUNDIALS library:


1. Configure with cmake, downloading and linking to SUNDIALS and NetCDF4:

   .. code-block:: bash

      cmake . -B build -DBOUT_DOWNLOAD_SUNDIALS=ON -DBOUT_DOWNLOAD_NETCDF_CXX4=ON

2. Build, compiling Hermes-3 and all dependencies using 4 parallel cores (adjust as necessary):

   .. code-block:: bash

      cmake --build build -j 4

3. Run the unit and integrated tests to check that everything is working:

   .. code-block:: bash

      cd build
      ctest

Note that the integrated tests require MPI, and so may not run on the
head nodes of many computing clusters.


Installing with SUNDIALS and PETSc (beuler)
~~~~~~~~~~

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
~~~~~~~~~~
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

.. _sec-slope-limiter-settings:

Slope (flux) limiter settings
~~~~~~~~~~

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
~~~~~~~~~~
Please see the `relevant page <https://bout-dev.readthedocs.io/en/stable/user_docs/advanced_install.html#optimisation-and-run-time-checking>`_ 
in the BOUT++ documentation.


Custom versions of BOUT++
~~~~~~~~~~

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
~~~~~~~~~~

The CMake configuration can be customised: See the `BOUT++
documentation
<https://bout-dev.readthedocs.io/en/latest/user_docs/installing.html#cmake>`_
for examples of using `cmake` arguments, or edit the compile options
interactively before building:

.. code-block:: bash

   ccmake . -B build


Troubleshooting issues
~~~~~~~~~~

The first step to troubleshooting compilation issues should always to delete
build folder for a fresh compilation. This can resolve several types of issues.

There have also been several reported issues due to Conda (e.g. making 
BOUT++ pick up the Conda MPI installation instead of the module one). A 
workaround is to compile with the CMake flag `-DBOUT_IGNORE_CONDA_ENV=ON`.


Using Spack
-----------

In this section we describe how to build Hermes-3 using `spack <https://spack.io/>`_ to manage the
installation of standard packages to your local environment. By default, dependencies like NetCDF4,
PETSc and SUNDIALS will be installed by spack, but note that it's also possible to `use your own
versions of packages <https://spack.readthedocs.io/en/latest/packages_yaml.html#external-packages>`_
(either system-installed or locally-built).

These instructions were last tested using Ubuntu 22.04.1 and spack version 0.23.1.
The default environment configuration assumes you have gcc installed.

Install Spack
~~~~~~~~~~~~~

Instructions for installing spack on a variety of operating systems can be found in the `spack docs <https://spack.readthedocs.io/en/latest/getting_started.html#installation>`_.
The commands below should work for most Debian-based Linux distributions.
First, install some prerequisites, e.g for Ubuntu:

.. code-block:: bash

   sudo apt update
   sudo apt install -y build-essential ca-certificates coreutils curl environment-modules gfortran git gpg lsb-release python3 python3-distutils python3-venv unzip zip

Now, clone spack (the command below checks out version 0.23.1).

.. code-block:: bash
   
   git clone -c feature.manyFiles=true --depth=2 https://github.com/spack/spack.git -b releases/v0.23

Initialise spack. Add this command to your .bashrc or similar to make spack available in all new shells:

.. code-block:: bash
  
   . spack/share/spack/setup-env.sh

The ``spack`` command should now be available; e.g.

.. code-block:: bash

   spack --version
    > 0.23.1 (2bfcc69fa870d3c6919be87593f22647981b648a)

Install Dependencies
~~~~~~~~~~~~~~~~~~~~

First, ensure you have the correct version of the `BOUT-spack` git submodule checked out:

.. code-block:: bash

   git submodule update --init BOUT-spack

At this point, standard spack commands can be used to activate the environment described in
``spack.yaml``. Alternatively, a wrapper script is provided that both activates the environment and
provides some useful bash functions and aliases. To use it, source the file with:

.. code-block:: bash

   . activate_h3env

You should see your prompt change to `[hermes-3]`, indicating that the spack environment is active.

.. note::

   The wrapper script runs ``spacktivate . -p -v gcc`` to load a 'view' when activating the
   environment. If you choose not to use the wrapper, you'll need to run a similar command in order
   for the instructions below to work.

.. note:: 

   If you've run ``spack install`` in this environment before, it's advisable to run ``spack
   concretize -f`` at this point to ensure the concretized 'spec' is up to date. See the `spack docs
   <https://spack.readthedocs.io/en/latest/environments.html>`_ for more details.

To install all dependencies:

.. code-block:: bash

   spack install --only dependencies -j 8

where the ``-j`` argument controls the number of parallel processes used to build packages.

This initial install takes some time to complete, because spack builds a large number of low-level
packages. It's possible to speed things up by defining a `packages.yaml
<https://spack.readthedocs.io/en/latest/packages_yaml.html>`_  that points to 'external' (system)
package versions, but unless storage space is a big concern, letting spack build its own versions is
usually the most trouble-free approach. This step rarely need to be repeated in its entirety unless
moving to another version of the same compiler, or switching to a different version of spack itself.

.. tip::
   Note that it's also possible to install hermes-3 itself as a package (remove ``--only
   dependencies`` from the above). This is the easiest approach if you just want to run the code,
   but re-installing it each time as part of a development workflow tends to be very slow. See below
   for an alternative.

Build Hermes-3 in the Spack environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A convenient way to build Hermes-3 in development is to use the dependencies installed via the
instructions above and then call ``cmake`` directly to configure and compile Hermes-3. The wrapper
script provides a bash function ``in_h3env`` that runs commands in the hermes-3 *build* environment,
setting all of the necessary paths to find headers, link libraries etc. To build
hermes-3 with CMake, run (e.g.):

.. code-block:: bash

   export h3_build="./builds/my_build" && in_h3env cmake -DBOUT_USE_PETSC=ON -B "$h3_build" && in_h3env cmake --build "$h3_build" -j8

In the above example, PETSc support is turned on at the configuration stage, which emulates the
result of installing hermes-3 directly with the default ``spec`` listed in spack.yaml
(``hermes-3%gcc+petsc``).

In order to run the tests, load the cmake package to make ``ctest`` available, then run as normal:

.. code-block:: bash

   spack load cmake
   cd ./builds/my_build
   ctest -j 3

If hermes-3 was installed with the +xhermes variant (as it is by default), the tests should pass
without having to modify any paths.

Configuration options
~~~~~~~~~~~~~~~~~~~~~

To see which `variants` of hermes-3 are available, run

.. code-block:: bash

   spack info --no-dependencies --no-versions hermes-3

By default, the top-level *spec* in the environment is `hermes-3%gcc+petsc`, which tells spack to
configure hermes-3 with PETSc support and build with gcc. To change this, first modify spack.yaml
(e.g.):

.. code-block:: yaml

   spack:
      specs:
         - hermes-3%gcc~petsc~xhermes
   ...

then reconcretize, and (re-)install packages as necessary:

.. code-block:: bash

   spack concretize -f
   spack install --only dependencies -j 8

.. warning:: 

   Note that the ``-f`` flag is crucial here. Without it, spack won't overwrite the existing spec
   (stored in spack.lock) and nothing will be updated.

.. important::
   If you're installing hermes-3 with CMake (i.e. using ``spack install --only dependencies``), this procedure will update the installed 
   dependencies according to the spec in spack.yaml, but **you'll still need to supply the appropriate configuration options to CMake**.
   If you're installing hermes-3 with spack directly, the correct configuration options will be set automatically.

Useful spack commands
~~~~~~~~~~~~~~~~~~~~~

See all dependencies required by the current spec:

.. code-block:: bash
  
   spack spec

Check paths to installation for lib, bin, and include files:

.. code-block:: bash
  
   spack find --paths [package_name]


Next Steps
----------

You are now ready to try running the examples in the ``build/examples/`` folder. See https://hermes3.readthedocs.io/en/latest/examples.html.