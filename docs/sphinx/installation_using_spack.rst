.. _sec-installation_using_spack:

Installation using Spack
========================

In these docs we describe how to install Hermes-3 with the assistance of Spack
to manage the installation of standard modules to the local environment.
More complicated modules like PETSc, SUNDIALS,
and BOUT++ that require configuration are installed by Hermes-3 automatically.
PETSc and SUNDIALS are available through
Spack, but further effort is required to understand
how to configure them correctly for Hermes-3
using the Spack interface. Theses installation instructions were tested on a
fresh Ubuntu 22.04 LTS.

Install Spack
-------------

See the ``spack`` docs here https://spack.readthedocs.io/en/latest/getting_started.html#installation.
First, get the required basic modules on your linux distribution. According to the Spack docs, these are as follows.

.. code-block:: bash

   sudo apt-get update
   sudo apt-get install build-essential ca-certificates coreutils curl environment-modules gfortran git gpg lsb-release python3 python3-distutils python3-venv unzip zip

Now, clone ``spack``.

.. code-block:: bash
   
   git clone -c feature.manyFiles=true --depth=2 https://github.com/spack/spack.git

Get the Spack functions onto the command line

.. code-block:: bash
  
   . spack/share/spack/setup-env.sh

Now you can use the ``spack`` function to manage your local environment. See https://spack.readthedocs.io/en/latest/basic_usage.html#installing-and-uninstalling.

Install required modules
------------------------

Install the required modules

.. code-block:: bash
   
   spack install cmake
   spack install fftw
   spack install openmpi
   spack install netcdf-c
   spack install netcdf-cxx4
   spack install python

Then load the modules (you may require to supply specific version numbers and hashes)

.. code-block:: bash
   
   spack load cmake
   spack load fftw
   spack load openmpi
   spack load netcdf-c
   spack load netcdf-cxx4
   spack load python

Check which modules are loaded with

.. code-block:: bash
  
   spack find --loaded


In a recent successful installation, the following modules were loaded.

.. code-block:: bash
  
   spack find --loaded
   -- linux-ubuntu22.04-skylake / gcc@11.4.0 -----------------------
   cmake@3.30.5  fftw@3.3.10  netcdf-c@4.9.2  netcdf-cxx4@4.3.1  openmpi@5.0.5  python@3.13.0
   ==> 6 loaded packages

Check paths to installation for lib, bin, and include files with

.. code-block:: bash
  
   spack find --paths module-of-interest

Make a virtual python environment with 

.. code-block:: bash
  
   python3 -m venv your-python-venv
   source /path/to/your-python-env/bin/activate

To install the required python libraries later.
You should install ``xbout`` https://github.com/boutproject/xBOUT 
and ``xhermes`` https://github.com/boutproject/xhermes with

.. code-block:: bash
  
   git clone https://github.com/boutproject/xBOUT.git
   cd xBOUT
   python3 -m pip install -e .
   git clone https://github.com/boutproject/xhermes.git
   cd xhermes
   python3 -m pip install -e .

Install PETSc
-------------

Download the latest PETSc, and configure it.

.. code-block:: bash
  
   wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.22.1.tar.gz
   tar -xf petsc-3.22.1.tar.gz
   cd petsc-3.22.1/ 
   ./configure --with-mpi=yes --download-hypre --download-make --with-fortran-bindings=0 --with-debugging=0 --download-fblaslapack=1

PETSc `configure` will now prompt you to make a command like

.. code-block:: bash
  
   make PETSC_DIR=/path/to/petsc-3.22.1/petsc-3.22.1 PETSC_ARCH=your-arch all

Check the install with

.. code-block:: bash
  
   make PETSC_DIR=/path/to/petsc-3.22.1/petsc-3.22.1 PETSC_ARCH=your-arch check

Export the appropriate environment variables

.. code-block:: bash
   
   export PETSC_DIR=/path/to/petsc-3.22.1/petsc-3.22.1
   export PETSC_ARCH=your-arch

Install Hermes-3
----------------

Now we are ready to install Hermes-3. First use 

.. code-block:: bash
  
   git clone https://github.com/bendudson/hermes-3.git
   cd hermes-3

Now run the configuration command

.. code-block:: bash
  
   cmake . -B build -DBOUT_DOWNLOAD_SUNDIALS=ON -DBOUT_USE_PETSC=ON

You will then be prompted to run

.. code-block:: bash
  
   cmake --build /home/mrhardman/hermes-3-work/hermes-3-spack/build

Test the install by 

.. code-block:: bash
   
   cd build
   ctest

To build in debug mode, use the flag ``-DCMAKE_BUILD_TYPE=Debug`` at the 
configuration step above.
Export a line like the following to your python path to make sure that 
python functions are available for post processing

.. code-block:: bash
 
   export PYTHONPATH=/path/to/hermes-3/build/external/BOUT-dev/tools/pylib:/path/to/hermes-3/external/BOUT-dev/tools/pylib:$PYTHONPATH

You are now ready to try running the example runs in the ``build/examples/`` folder. See https://hermes3.readthedocs.io/en/latest/examples.html.
