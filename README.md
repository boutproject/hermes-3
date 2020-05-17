      _   _                                   ______
     | | | |                                 |___   |
     | |_| | ___ _ __ _ __ ___   ___  ___        / /
     |  _  |/ _ \ '__| '_ ` _ \ / _ \/ __|      / /
     | | | |  __/ |  | | | | | |  __/\__ \     / /___
     \_| |_/\___|_|  |_| |_| |_|\___||___/    |______| 


Hermes plasma edge simulation model. Uses BOUT++ framework, adds finite volume
operators and neutral gas models.

This is Hermes-3, a hot ion multifluid drift-reduced model.

Author: Ben Dudson, University of York <benjamin.dudson@york.ac.uk>

Released under the GPL license

License
=======

Full text of the license is in the file LICENSE. If you are using Hermes-2,
please cite the relevant papers.

    Copyright B.Dudson, J.Leddy, University of York, September 2017-2020
              email: benjamin.dudson@york.ac.uk

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

Installing BOUT++
=================

This version works with BOUT++ v4.3 or later

    git clone https://github.com/boutproject/BOUT-dev.git
    cd BOUT-dev

To run this model, preconditioning is strongly recommended, and requires the CVODE solver, part of [SUNDIALS](http://computation.llnl.gov/projects/sundials).
Tested with version 2.6.0. To enable CVODE, BOUT++ should be configured using

    ./configure --with-cvode

or

    ./configure --with-sundials

(which then also enables the IDA solver). Compile BOUT++ with

    make

Compiling Hermes
================

To compile, run "make" and specify the location of BOUT++
> $ make BOUT_TOP=/path/to/BOUT/

Testing
=======

> $ make check BOUT_TOP=/path/to/BOUT

This will run both unit and integrated tests.

Examples
========

> $ ./hermes-3 -d examples/blob2d

