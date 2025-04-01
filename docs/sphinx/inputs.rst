.. _sec-configuration:

Input and grid files
===========

Inputs
-----------

Hermes-3 simulations are configured using a BOUT.inp "options"
file. This file has sections marked with square brackets
(e.g. ``[section]``), and ``key = value`` options within sections. The
format is described in the `BOUT++ options manual
<https://bout-dev.readthedocs.io/en/stable/user_docs/bout_options.html>`_.

The structure of the input file is divided into three sections.
The first section contains top-level BOUT++ settings such as the 
number of output timesteps ``nout`` and the output timestep size
``timestep`` (in normalised units where 95788 is one milisecond).
Note that the solver timestep is adaptive and not user-settable.

This is followed by ``[mesh]``, ``[solver]`` and ``[hermes]`` headers, where the ``[hermes]``
section defines the list of components used. The component order
matters, as the components are executed in order.

.. code-block:: ini

    nout = 200
    timestep = 5000

    [mesh]
    ...

    [solver]
    ...

    [hermes]
    components = (...)

Top-level components are defined underneath the ``[hermes]`` section:

.. code-block:: ini

    [hermes]
    components = (...)

    [sheath_boundary_simple]
    gamma_i = 3.5  # Ion sheath heat transmission coefficient
    gamma_e = 4.8  # Electron sheath heat transmission coefficient

    [sound_speed]
    electron_dynamics = false

    [d+]
    type = (evolve_density, evolve_pressure, evolve_momentum,
        noflow_boundary, upstream_density_feedback)

A top-level component is either a component which requires to have information
regarding all species (e.g. ``[sheath_boundary_simple]`` or ``[sound_speed]``), 
or a species itself (e.g. ``[d+]``). Apart from species, each component name
has a corresponding header and/or implementation file, e.g.
``sheath_boundary_simple.hxx`` and ``sheath_boundary_simple.cxx``.

Species-level components are instantiated individually for each species, and
so are defined inside the top-level species component. Common species-level components
are ``evolve_density``, ``evolve_momentum`` and ``evolve_pressure`` which correspond
to the continuity, momentum and pressure equations, respectively. This is also 
where other top-level components such as ``recycling`` have their options set.

This is followed by sections for individual evolved variables, e.g. ``Nd+``,
``Pd+`` and ``NVd+``. These allow for the setting of the initial conditions (``function``),
sources (``source``) and boundary conditions (e.g. ``bndry_all``).

The conventions and syntax for defining the initial conditions and source are
covered in the `BOUT++ documentation 
<https://bout-dev.readthedocs.io/en/stable/user_docs/variable_init.html>`_`.
Simple boundary conditions are defined in the `boundary conditions 
section <https://bout-dev.readthedocs.io/en/stable/user_docs/boundary_options.html>`_.



Gridding
--------------
Hermes-3 can be configured using an analytically defined grid, or one created
in Hypnotoad.


Guard cells
~~~~~~~~~~~~~~


Grid redistribution and interpolation
~~~~~~~~~~~~~~


Using grids to pass inputs to Hermes-3
~~~~~~~~~~~~~~


Metric coefficients
~~~~~~~~~~~~~~

The option ``hermes:recalculate_metric`` controls how the metric tensor is calculated. 
By default ``recalculate_metric`` is ``false``, meaning that the metric tensor
components (``g11``, ``g_22`` etc.) are taken from the grid file.

Setting ``recalculate_metric`` to ``true`` causes Hermes-3 to read
``Rxy``, ``Bpxy`` and other geometric quantities from the grid file.
The metric tensor is recalculated to the orthogonal field-aligned
coordinate system described in the `BOUT++ coordinate manual
<https://bout-dev.readthedocs.io/en/stable/user_docs/coordinates.html#jacobian-and-metric-tensors>`_.

**Note** Previous Hermes-3 versions had an option ``loadmetric`` with
the same behavior but the opposite default (``loadmetric=false``
rather than ``recalculate_metric=true``).


If ``hermes:recalculate_metric`` is false (the default), then the coordinate
metrics loaded from the grid file are usually in SI units.  By default
``normalise_metric`` is ``true``, and the loaded metrics are
normalised using the Hermes-3 normalisation factors.

If ``recalculate_metric`` is set to ``true`` then the metrics will always
be normalised, and the ``normalise_metric`` option is not used.
The default BOUT++ behavior is to throw an exception if an option is
set but not used.

