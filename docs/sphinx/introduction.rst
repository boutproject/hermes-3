.. _sec-introduction:

Introduction
============

Hermes-3 is a plasma simulation model built on `BOUT++
<http://boutproject.github.io/>`_, developed mainly for simulating the
edge of magnetically confined plasmas such as tokamaks. The source
code is `available on Github
<https://github.com/boutproject/hermes-3>`_. The main aim of this model
is multi-species simulation of fusion reactors, where the plasma will
contain a mixture of deuterium, tritium, helium and other species.

Hermes-3 is multi-fidelity, allowing the simulation of 1D, 2D and 3D 
tokamak plasmas both in steady-state and unsteady / turbulent regimes.

It is designed to be modular and organised into reusable
components, which can be tested individually and then configured at
run-time.

NOTE:
  This is a research code and may change without notice. The documentation
  is a work in progress and may not be complete. Please raise any issues
  on the `GitHub repo <https://github.com/boutproject/hermes-3>`_.


Publications
-------------

* B.Dudson, M.Kryjak, H.Muhammed, P.Hill, J,Omotani `Hermes-3:
  Multi-component plasma simulations with
  BOUT++ <https://doi.org/10.1016/j.cpc.2023.108991>`_
  Comp. Phys. Comm. 2023
  108991. doi: `10.1016/j.cpc.2023.108991 <https://doi.org/10.1016/j.cpc.2023.108991>`_.
  Preprint:
  `arXiv.2303.12131 <https://doi.org/10.48550/arXiv.2303.12131>`_.

* G.K. Holt, A. Keats, S. Pamela, M. Kryjak, A. Agnello,
  N.C. Amorisco, B.D. Dudson and M. Smyrnakis `Tokamak divertor plasma
  emulation with machine
  learning <https://doi.org/10.1088/1741-4326/ad4f9e>`_ 2024
  Nucl. Fusion 64 086009
  doi: `10.1088/1741-4326/ad4f9e <https://doi.org/10.1088/1741-4326/ad4f9e>`_

* Thomas Body, Thomas Eich, Adam Kuang, Tom Looby, Mike Kryjak, Ben Dudson, Matthew Reinke
  `Detachment scalings derived from 1D scrape-off-layer simulations <https://doi.org/10.1016/j.nme.2024.101819>`_
  Nucl. Mat. Energy 2024 101819
  doi: `10.1016/j.nme.2024.101819 <https://doi.org/10.1016/j.nme.2024.101819>`_

* Huayi Chang, Ben Dudson, Jizhong Sun, Mike Kryjak, Yang Ye, Mao Li,
  Weikang Wang `Hermes-3 simulation of the low-n X-point mode driven
  by impurity in tokamak edge
  plasmas <https://doi.org/10.1016/j.nme.2025.101913>`_ Nucl. Mat. Energy 2025 101913
  doi: `10.1016/j.nme.2025.101913 <https://doi.org/10.1016/j.nme.2025.101913>`_


FAQs
-------------

How do I join the Hermes-3 community?
  We have regular group meetings and a Slack channel.
  Contact Mike (mike.kryjak@ukaea.uk) or Ben (dudson2@llnl.gov) for details.

How do I install Hermes-3?
  Should be relatively straightforward, whether using CMake or Spack.
  If compiling without a module environment (e.g. on a laptop),
  Spack is recommended. See :ref:`sec-installation` for details.

How do I configure a simulation?
  Hermes-3 uses input files formatted in a Python-like syntax,
  see :ref:`sec-inputs`. Use the provided 
  :ref:`sec-examples` cases as a starting point. Check 
  the ``hermes-3\examples\`` directory for extra examples.
  Recommended solver settings are in the :ref:`sec-solver_numerics` section.

How do I run Hermes-3, how do I restart simulations?
  See :ref:`sec-execution`.

How do I post-process a simulation?
  There are two interfaces: the legacy Python package boutdata
  and the newer, Xarray powered xHermes. xHermes can be a bit 
  slower but handles normalisation for you and has many quality
  of life features. See section :ref:`sec-postprocessing`.

What equations am I solving?
  Hermes-3 is modular and the solved equations are built from "components".
  Please refer to sections :ref:`sec-equations` and :ref:`sec-closure`.
  There is a separate section on :ref:`sec-reactions`, including impurity
  radiation, and one on :ref:`sec-boundary_conditions`. 

What about tests?
  See :ref:`sec-tests` for details on the test suite.

How can I modify Hermes-3?
  Please see the :ref:`sec-developer` section for details on how to
  contribute to Hermes-3. We welcome contributions, whether they are
  bug fixes, new features or documentation improvements.


License
------------

Hermes-3 is released under the GPL-3 license. If you are using Hermes-3, please
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



