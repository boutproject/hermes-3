# Hermes-3

[![License](https://img.shields.io/badge/license-GPL-blue.svg)](https://img.shields.io/badge/license-GPL-blue.svg)
![Build status](https://github.com/bendudson/hermes-3/workflows/Tests/badge.svg)

Hermes-3 is a multifluid plasma simulation model for transport and
turbulence in the edge of magnetically confined plasmas, such as
tokamaks. It is built on the [BOUT++
framework](https://github.com/boutproject/BOUT-dev/), and uses a
system of reusable components to build models at runtime based on
input configuration, in 1D, 2D or 3D curvlinear coordinates.  The
manual is [here on
Readthedocs](https://hermes3.readthedocs.io/en/latest/).

Post-processing Hermes-3 is made easy by using the [xHermes](https://github.com/boutproject/xhermes) Python library.
xHermes is built on top of [xBOUT](https://github.com/boutproject/xBOUT) and Xarray and automatically accounts 
for the normalisation of Hermes-3 parameters, as well as providing access to many useful xBOUT post-processing functions.
Many of these are showcased in the [xBOUT examples](https://github.com/boutproject/xBOUT-examples).

*Note* Under development, research code, may change without notice.

## License

Hermes-3 is released under the GPL-3 license. See [LICENSE](./LICENSE)
and [NOTICE](./NOTICE) for details.  If you are using Hermes-3, please
cite the relevant papers.

All new contributions must be made under the GPLv3 license.

LLNL-CODE-845139


    Copyright Hermes-3 contributors 2017-2025
              email: dudson2@llnl.gov

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


## Installing and testing

Only CMake is supported for building Hermes-3 and running the tests.
During configuration
[BOUT++](https://github.com/boutproject/BOUT-dev/) will be
automatically downloaded as a submodule, together with some
dependencies. [NetCDF](https://www.unidata.ucar.edu/software/netcdf/)
and [FFTW](https://www.fftw.org/) are assumed to be installed already;
optional dependencies include
[SUNDIALS](https://computing.llnl.gov/projects/sundials) and
[PETSc](https://petsc.org). The recommended way to build Hermes-3
links to the [SUNDIALS](https://computing.llnl.gov/projects/sundials)
library.

1) Configure with cmake, downloading and linking to SUNDIALS:

    $ cmake . -B build -DBOUT_DOWNLOAD_SUNDIALS=ON

2) Build, compiling Hermes-3 and all dependencies:

    $ cmake --build build

3) Run the unit and integrated tests to check that everything is working:

    $ cd build
    $ ctest

Note that the integrated tests require MPI, and so may not run on the
head nodes of many computing clusters.

The CMake configuration can be customised: See the [BOUT++
documentation](https://bout-dev.readthedocs.io/en/latest/user_docs/installing.html#cmake)
for examples of using `cmake` arguments, or edit the compile options
interactively before building:

    $ ccmake . -B build

If you have already installed BOUT++ and want to use that rather than
configure and build BOUT++ again, set `HERMES_BUILD_BOUT` to `OFF` and pass
CMake the path to the BOUT++ `build` directory e.g.

    $ cmake . -B build -DHERMES_BUILD_BOUT=OFF -DCMAKE_PREFIX_PATH=$HOME/BOUT-dev/build

Note that Hermes-3 currently requires a specific version of BOUT++:
https://github.com/boutproject/BOUT-dev/commit/7152948fbde505f6708d5ca4a9c21e5828d1e0a1

## Examples

There are example inputs under the examples/ subdirectory. A simple
example is a 2D (drift plane) simulation of a plasma blob/filament,
similar to the BOUT++
[blob2d](https://github.com/boutproject/BOUT-dev/tree/master/examples/blob2d)
example:

    ./hermes-3 -d examples/blob2d

See the
[examples](https://github.com/bendudson/hermes-3/tree/master/examples)
for more complicated cases.

## Publications

* B.Dudson, M.Kryjak, H.Muhammed, P.Hill, J,Omotani [Hermes-3:
  Multi-component plasma simulations with
  BOUT++](https://doi.org/10.1016/j.cpc.2023.108991)
  Comp. Phys. Comm. 2023
  108991. doi:[10.1016/j.cpc.2023.108991](https://doi.org/10.1016/j.cpc.2023.108991).
  Preprint:
  [arXiv.2303.12131](https://doi.org/10.48550/arXiv.2303.12131).
* Thomas Body, Thomas Eich, Adam Kuang, Tom Looby, Mike Kryjak, Ben Dudson, Matthew Reinke
  [Detachment scalings derived from 1D scrape-off-layer simulations](https://doi.org/10.1016/j.nme.2024.101819)
  Nucl. Mat. Energy 2024 101819
  doi:[10.1016/j.nme.2024.101819](https://doi.org/10.1016/j.nme.2024.101819)
* Huayi Chang, Ben Dudson, Jizhong Sun, Mike Kryjak, Yang Ye, Mao Li,
  Weikang Wang [Hermes-3 simulation of the low-n X-point mode driven
  by impurity in tokamak edge
  plasmas](https://doi.org/10.1016/j.nme.2025.101913) Nucl. Mat. Energy 2025 101913
  doi:[10.1016/j.nme.2025.101913](https://doi.org/10.1016/j.nme.2025.101913)
