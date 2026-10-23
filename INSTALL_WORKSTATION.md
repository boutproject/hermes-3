# Installing Hermes-3 on an Ubuntu workstation

How to build Hermes-3 with PETSc on Ubuntu 24.04, using the system OpenMPI and netCDF.
This setup is known to work on `titipounamu`; the versions used there are noted at
each step.

The build uses **system libraries only**. The conda/mamba `hermes` environment is for
Python post-processing and tests; keep it out of the build (step 4 shows how, even with
the environment active).

## 1. System packages

```bash
sudo apt install build-essential gfortran git \
    libopenmpi-dev openmpi-bin \
    libnetcdf-dev libnetcdf-c++4-dev libhdf5-dev \
    libfftw3-dev libopenblas-dev
sudo snap install cmake --classic      # needs CMake >= 3.23
```

Versions on `titipounamu`: GCC 13.3, Open MPI 4.1.6, netCDF-C 4.9.2, netCDF-C++ 4.3.1,
HDF5 1.10.10 (serial), CMake 4.4.

Check that `mpicc`/`mpicxx` are Open MPI, not MPICH:

```bash
mpicxx --showme:version      # should print "Open MPI 4.1.x"
```

## 2. PETSc (built against Open MPI)

PETSc must use the **same MPI as Hermes-3** (the system Open MPI). A PETSc built with
MPICH compiles but fails to link or run.

```bash
cd ~
git clone -b release https://gitlab.com/petsc/petsc.git petsc
cd petsc
./configure --prefix=$HOME/petsc-build --with-mpi=yes \
    --with-debugging=0 --with-fortran-bindings=0 --with-make-np=32 \
    --download-hypre --download-mumps --download-strumpack --download-scalapack \
    --download-superlu --download-superlu_dist --download-ptscotch \
    --download-parmetis --download-metis --download-zfp \
    --download-slepc --download-hpddm --download-bison --download-openblas=1
make all          # use the exact command configure prints at the end
make install      # copies everything, including the downloaded libraries, to $HOME/petsc-build
```

Version on `titipounamu`: PETSc 3.26.0.

Use the **installed** prefix (`$HOME/petsc-build`) in step 4, not the in-tree
`petsc/arch-*` directory: the downloaded packages (HYPRE, MUMPS, SCOTCH, ...) are only
complete in the installed prefix. Check that nothing is missing:

```bash
ldd ~/petsc-build/lib/libpetsc.so | grep 'not found'    # must print nothing
ldd ~/petsc-build/lib/libpetsc.so | grep libmpi.so      # must be the system Open MPI
```

## 3. Hermes-3 source

```bash
cd ~
git clone --recursive https://github.com/boutproject/hermes-3.git
cd hermes-3
git checkout dev-OS-dipole          # the dipole branch
git submodule update --init --recursive
```

The submodules (`external/BOUT-dev`, `external/BOUT-spack`, `external/googletest`) are
required: BOUT++ is built as part of Hermes-3.

## 4. Configure and build

```bash
cd ~/hermes-3
cmake . -B build \
    -DBOUT_DOWNLOAD_SUNDIALS=ON \
    -DBOUT_DOWNLOAD_NETCDF_CXX4=OFF \
    -DBOUT_USE_PETSC=ON \
    -DPETSC_DIR=$HOME/petsc-build \
    -DHERMES_SLOPE_LIMITER=VanAlbada \
    -DnetCDF_DIR=/usr/lib/x86_64-linux-gnu/cmake/netCDF \
    -UHDF5_DIR \
    -DCMAKE_IGNORE_PREFIX_PATH=$HOME/.conda/envs/hermes
cmake --build build -j32
```

What each of the non-obvious options is for:

| Option | Why |
|---|---|
| `-DnetCDF_DIR=/usr/...` | Pins netCDF-C to the system copy, which matches the system netCDF-C++ |
| `-UHDF5_DIR` | Clears a conda HDF5 path cached by an earlier configure |
| `-DCMAKE_IGNORE_PREFIX_PATH=...conda...` | Lets you configure with the `hermes` env active without CMake picking up its netCDF, HDF5 or curl |
| `-DHERMES_SLOPE_LIMITER=VanAlbada` | Smooth (differentiable) advection limiter; converges better with the Newton (SNES) solver than the default `MC` |

Check that only system netCDF libraries are linked (no `.conda` path):

```bash
tr ' ' '\n' < build/CMakeFiles/hermes-3.dir/link.txt | grep -E 'netcdf|conda'
```

The executable is `build/hermes-3`. On `titipounamu`, `hermes-3-dev` in
`/usr/local/bin` is a symlink to it, so every rebuild changes what `hermes-3-dev` runs.
Jobs that are already running keep their old binary.

## 5. Python environment (post-processing and tests)

```bash
mamba create -n hermes -c conda-forge python numpy scipy matplotlib xarray netcdf4
mamba activate hermes
pip install boutdata boututils xbout
```

Versions on `titipounamu`: Python 3.14, numpy 2.4, xarray 2026.2, netCDF4 1.7.4,
boutdata/boututils 0.3.0, xbout 0.3.8.

## 6. Test

```bash
cd ~/hermes-3/build
./hermes_unit_tests --gtest_brief=1      # expect: all tests PASSED
ctest -j8 --output-on-failure            # run with the hermes env active (needs Python)
```

Known `ctest` failures with this Python environment (Python 3.14, recent xarray), all in
the Python post-processing, not in Hermes-3 itself: `2D-production`, `2D-recycling`,
`1D-recycling`, `1D-recycling-dthe`. `examples` fails if a previous interrupted run left
`build/examples/tokamak-2D/json_database` behind; delete it and rerun.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `libhdf5_serial.so.103: undefined reference to 'curl_...@CURL_OPENSSL_4'` | Configured with the conda env active: conda netCDF mixed with system netCDF-C++, and conda's `libcurl` found at link time | Re-run step 4 exactly (it sets `netCDF_DIR`, clears `HDF5_DIR`, ignores the conda prefix) |
| `fatal error: petsc.h: No such file or directory` | The PETSc install that `build/` was configured with was moved or deleted | Reinstall PETSc (step 2) and re-run step 4 |
| `libpetsc.so: undefined reference to 'HYPRE_...'` or `'SCOTCH_...'` | `PETSC_DIR` points at the in-tree `petsc/arch-*` build | Use the installed prefix (`$HOME/petsc-build`) |
| Link or runtime MPI errors after a PETSc rebuild | PETSc built with MPICH | Rebuild PETSc with the system Open MPI (check `MPICC_SHOW` in `~/petsc-build/lib/petsc/conf/petscvariables`) |
| `hermes-3-dev: No such file or directory` after a failed build | A failed link deletes the old `build/hermes-3` before writing the new one | Fix the build error and rebuild |
