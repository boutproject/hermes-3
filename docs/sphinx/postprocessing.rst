.. _sec-postprocessing:

Post-processing
-----------

Tools
~~~~~~~~~~~

There are two workflows for post-processing Hermes-3 results. The legacy way
is to use the `boutdata <https://github.com/boutproject/boutdata>`_ package, which
can be used to extract the raw variables from the dump files individually, which
then have to be unnormalised to be converted into SI units.

The recommended way is to use `xHermes <https://github.com/boutproject/xhermes>`_, 
which is built upon the Xarray BOUT++ post-processing tool `xBOUT 
<https://github.com/boutproject/xBOUT>`_.
The main advantage of xHermes is that it contains several routines for pre-processing
the data, including automatic unnormalisation and geometry calaculations. The disadvantage
is that it may be slower than ``boutdata`` when running in parallel due to the netCDF
I/O overhead. This can be mitigated by squashing the dump files. Efforts to improve
the performance are ongoing.

Please refer to the xHermes `readme <https://github.com/boutproject/xhermes?tab=readme-ov-file#xhermes>`_
and `examples <https://github.com/boutproject/xhermes/tree/main/examples>`_ for details.

List of diagnostic variables
~~~~~~~~~~~

It can be problematic to record all of the diagnostic variables in the documentation
as the code is still in active development. Thankfully, all of the diagnostics are
saved with helpful metadata containing the name, units, normalisation, description
and the origin of the variable from within the code.

This is available in xHermes. You can access it by the following example:

.. code-block:: ini

   ds["Nd"].attrs()

The attributes are:

time_dimension
   This indicates that the variable is time evolving. If you remove this line,
   it will be saved only at the first RHS evaluation.

units
   A string showing the units for post-processing. xHermes picks this up.

conversion
   A float representing the normalisation factor. xHermes picks this up to do
   automatic conversion to SI units.

standard_name, long_name
   We aren't consistent on what should be in each, but they are meant to describe
   the variables in post-processing.

species, source
   The relevant species and component that the diagnostic is coming from

Alternatively, the same attributes can be seen in the source code.


Reaction diagnostics convention
~~~~~~~~~~~

There are five categories of diagnostics. ``K`` is a simple reaction rate (always positive).
``S`` is a source of density and is positive if the ion density is increasing.
``E`` is a source of energy due to energy transfer (e.g. charge exchange) and 
is positive if the ion energy is increasing.
``F`` is a momentum transfer rate and is currently negative if the ion energy is increasing.
``R`` is a loss of energy from the system, e.g. due to radiation. It is always negative.

+------------------+---------------------------+-------------------------+
| Variable prefix  |   Units                   | Description             |
+==================+===========================+=========================+
| K                |   :math:`s^{-1}`          | Reaction rate           |
+------------------+---------------------------+-------------------------+
| S                |   :math:`m^{-3}s^{-1}`    | Density transfer source |
+------------------+---------------------------+-------------------------+
| E                |   :math:`Wm^{-3}`         | Energy transfer source  |
+------------------+---------------------------+-------------------------+
| F                |   :math:`kgm^{-2}s^{-2}`  | Momentum transfer rate  |
+------------------+---------------------------+-------------------------+
| R                |   :math:`Wm^{-3}`         | System loss source      |
+------------------+---------------------------+-------------------------+

Charge exchange can be slightly complicated:

For two species `a` and `b`, the channel `Fab_cx` is a source of momentum for species `a` due to
charge exchange with species `b`. There are corresponding sinks for
the products of the charge exchange reaction which are not saved.

For example,reaction `d + t+ -> d+ + t` will save the following
forces (momentum sources):
- `Fdt+_cx` is a source of momentum for deuterium atoms `d` and sink of momentum for deuterium ions `d+`.
- `Ft+d_cx` is a source of momentum for tritium ions `t+` and sink of momentum for tritium atoms `t`.

The reason for this convention is the existence of the inverse reactions:
`t + d+ -> t+ + d` outputs diagnostics `Ftd+_cx` and `Fd+t_cx`.

Flow diagnostics convention
~~~~~~~~~~~

Hermes-3 saves a number of different flow rates at a cell boundary. The boundary is 
located at either the "xlow" or "ylow" side, where "low" refers to the negative 
direction in index space, and "up" would refer to the positive direction. For example,
the "xlow" side of the cell at (i,j,k) is the boundary between cells (i,j,k) and (i-1,j,k).
Positive values correspond to flow into the cell.

There are several types of flow diagnostics at the moment:

+------------------+---------------------------+-------------------------+
| Variable prefix  |   Units                   | Description             |
+==================+===========================+=========================+
| pf               |   :math:`s^{-1}`          | Particle flow           |
+------------------+---------------------------+-------------------------+
| ef               |   :math:`W`               | Energy flow             |
+------------------+---------------------------+-------------------------+
| mf               |   :math:`N`               | Momentum flow           |
+------------------+---------------------------+-------------------------+

There are numerous variants of the above throughout the code, which are
described in their metadata. Note that the diagnostics have not yet been 
implemented for each term.



