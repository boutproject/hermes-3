.. _sec-introduction:

Introduction
============

Hermes-3 is a plasma simulation model built on `BOUT++
<http://boutproject.github.io/>`_, developed mainly for simulating the
edge of magnetically confined plasmas such as tokamaks. The source
code is `available on Github
<https://github.com/bendudson/hermes-3>`_. The main aim of this model
is multi-species simulation of fusion reactors, where the plasma will
contain a mixture of deuterium, tritium, helium and other species.

The unique feature of Hermes-3 is that it is modular and organised into reusable
components, which can be tested individually and then configured at
run-time.

FAQs
-------------

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

How do I post-process a simulation?
  There are two interfaces: the legacy Python package boutdata
  and the newer, Xarray powered xHermes. xHermes can be a bit 
  slower but handles normalisation for you and has many quality
  of life features. See section :ref:`sec-execution-postprocessing`.

What equations am I solving?
  Hermes-3 is modular and the solved equations are built from "components".
  Please refer to sections :ref:`sec-equations` and :ref:`sec-closure`.
  There is a separate section on :ref:`sec-reactions`, including impurity
  radiation, and one on :ref:`sec-boundary_conditions`. 

What about tests?
  See :ref:`sec-developer` for details on the test suite.

How do I contribute?
  Please see the :ref:`sec-developer` section for details on how to
  contribute to Hermes-3. We welcome contributions, whether they are
  bug fixes, new features or documentation improvements.


