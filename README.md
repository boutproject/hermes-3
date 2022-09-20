# Hermes-3

[![License](https://img.shields.io/badge/license-GPL-blue.svg)](https://img.shields.io/badge/license-GPL-blue.svg)
![Build status](https://github.com/bendudson/hermes-3/workflows/Tests/badge.svg)

Hermes plasma edge simulation model. Uses BOUT++ framework, adds finite volume
operators and neutral gas models.

This is Hermes-3, a hot ion multifluid drift-reduced model. The manual is
[here on Readthedocs](https://hermes3.readthedocs.io/en/latest/).

*Note* Under development, research code, may change without notice.

Released under the GPL license

## License

Full text of the license is in the file LICENSE. If you are using Hermes-3,
please cite the relevant papers.

    Copyright Hermes-3 contributors 2017-2022
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

### CMake

To use the default options and minimal dependencies just run:

    $ cmake . -B build

Alternatively the CMake build can be customised: See the [BOUT++
documentation](https://bout-dev.readthedocs.io/en/latest/user_docs/installing.html#cmake)
for examples of using `cmake` arguments, or edit the compile options
interactively before building:

    $ ccmake . -B build

During configuration
[BOUT++](https://github.com/boutproject/BOUT-dev/) will be
automatically downloaded as a submodule, together with some
dependencies (NetCDF and FFTW are assumed to be installed already,
along with optional dependencies like SUNDIALS and PETSc if they are
requested).

If you have already installed BOUT++ and want to use that rather than
configure and build BOUT++ again, set `HERMES_BUILD_BOUT` to `OFF` and pass
CMake the path to the BOUT++ `build` directory e.g.

    $ cmake . -B build -DHERMES_BUILD_BOUT=OFF -DCMAKE_PREFIX_PATH=$HOME/BOUT-dev/build

Once configured, run build to compile BOUT++ and then
Hermes-3:

    $ cmake --build build

Then run the unit and integrated tests to check that everything is working:

    $ cd build
    $ ctest

### Autotools

*Note*: Building with `configure` & `make` is no longer supported, and
may be completely removed in future. The CMake build ensures that the
version of BOUT++ downloaded is compatible with this version of
Hermes, whereas the procedure below does not.

To build and run tests with GNU autoconf and make, first install BOUT++:

    git clone -b next https://github.com/boutproject/BOUT-dev.git BOUT-next
    cd BOUT-next

To run some cases, preconditioning is strongly recommended, and
requires the CVODE solver, part of
[SUNDIALS](http://computation.llnl.gov/projects/sundials).
To enable CVODE, BOUT++ should be configured using

    ./configure --with-cvode

or

    ./configure --with-sundials

(which then also enables the IDA solver). Compile BOUT++ with

    make


Then clone the Hermes-3 repository

    git clone https://github.com/bendudson/hermes-3

    cd hermes-3

To compile, run "make" and specify the location of the BOUT++
installation

    make BOUT_TOP=/path/to/BOUT-next

This path should be the full path, not relative path, to avoid
problems with compilation in subdirectories.

To run the tests, which are run on Travis:

    make check BOUT_TOP=/path/to/BOUT-next

This will run both unit and integrated tests.

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
