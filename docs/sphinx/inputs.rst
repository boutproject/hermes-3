.. _sec-inputs:

Inputs
======

Hermes-3 simulations are configured using a BOUT.inp "options"
file. This file has sections marked with square brackets
(e.g. ``[section]``), and ``key = value`` options within sections. The
format is described in the `BOUT++ options manual
<https://bout-dev.readthedocs.io/en/stable/user_docs/bout_options.html>`_.

Below are described some of the important options in the BOUT.inp
file for Hermes-3 simulations.

hermes section
--------------

The [hermes] section of the BOUT.inp file is the starting point for
configuring the system of equations that Hermes-3 will solve.


hermes:recalculate_metric
~~~~~~~~~~~~~~~~~~~~~~~~~

This option controls how the metric tensor is calculated. By default
``recalculate_metric`` is ``false``, meaning that the metric tensor
components (``g11``, ``g_22`` etc.) are taken from the grid file.

Setting ``recalculate_metric`` to ``true`` causes Hermes-3 to read
``Rxy``, ``Bpxy`` and other geometric quantities from the grid file.
The metric tensor is recalculated to the orthogonal field-aligned
coordinate system described in the `BOUT++ coordinate manual
<https://bout-dev.readthedocs.io/en/stable/user_docs/coordinates.html#jacobian-and-metric-tensors>`_.

**Note** Previous Hermes-3 versions had an option ``loadmetric`` with
the same behavior but the opposite default (``loadmetric=false``
rather than ``recalculate_metric=true``).

hermes:normalise_metric
~~~~~~~~~~~~~~~~~~~~~~~

If ``recalculate_metric`` is false (the default), then the coordinate
metrics loaded from the grid file are usually in SI units.  By default
``normalise_metric`` is ``true``, and the loaded metrics are
normalised using the Hermes-3 normalisation factors.

If ``recalculate_metric`` is set to ``true`` then the metrics will always
be normalised, and the ``normalise_metric`` option is not used.
The default BOUT++ behavior is to throw an exception if an option is
set but not used.
