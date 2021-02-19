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

An unusual feature of this model is that it is organised into reusable
components, which can be tested individually and then configured at
run-time. For example a transport simulation with deuterium and tritium ions and
atoms has an input file specifying the components

.. code-block:: ini
  
  [hermes]
  components = d+, d, t+, t, e, collisions, sheath_boundary, recycling, reactions

The governing equations for each species are specified e.g.

.. code-block:: ini

  [d+]
  type = evolve_density, evolve_momentum, evolve_pressure, anomalous_diffusion
  AA = 2   # Atomic mass
  charge = 1

and other components have their configuration options e.g. for reactions

.. code-block:: ini

  [reactions]
  type = (
          d + e -> d+ + 2e,   # Deuterium ionisation
          t + e -> t+ + 2e,   # Tritium ionisation
         )

