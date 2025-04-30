.. _sec-geometry_1d:

Fieldline Geometry Component
===========================

The ``fieldline_geometry`` component in Hermes-3 helps to set up 1D simulations by providing a convenient way to modify the effective area of the flux tube.

Functionality
-----------

This component calculates and applies broadening effects to the flux tube, effectively altering its cross-sectional area along the parallel direction. It considers three main factors:

* ``geometric_broadening``:  :math:`f_{geo}=R(L_\parallel)/R_u` accounts for the increase in the major radius along the flux tube. 
* ``transport_broadening``:  :math:`f_{transp}=\lambda_{INT}(L_\parallel)/\lambda_{q,u}` represents the cross-field transport in the divertor (i.e. the :math:`S` term in :math:`\lambda_{INT}`).
* ``flux_expansion``: :math:`f_{exp}=\sin(\theta)_u/\sin(\theta)(L_\parallel)` with :math:`\sin(\theta)=B_{pol}/B_{tot}` captures the variation in the magnetic field pitch angle.

These effects are combined to calculate a total broadening factor, which is then used to modify the Jacobian :math:`J=1/B`.

Implementation
--------------

The component is implemented in ``hermes-3/src/fieldline_geometry.cxx`` and its header file is ``hermes-3/include/fieldline_geometry.hxx``.

Usage
-----

To utilize this component, it must be included in the ``components`` list within the ``BOUT.inp`` file. For example:

.. code-block:: ini

 [hermes]
 components = (fieldline_geometry, d+, d, e,
 sheath_boundary_simple, collisions, recycling, reactions,
 electron_force_balance, neutral_parallel_diffusion)

The component is configured by adding a ``[fieldline_geometry]`` block to the ``BOUT.inp`` file.  In this block, the broadening functions can be defined.  For example:

.. code-block:: ini

 [fieldline_geometry]
 transport_broadening = where({lpar} > mesh:length_xpt, 2.0, 1.0)
 flux_expansion = 1.0
 diagnose = true

* If a broadening variable (``transport_broadening``, ``flux_expansion``, or ``geometric_broadening``) is not set, it defaults to 1 (constant).
* The parallel length from the stagnation point is accessed using the variable ``{lpar}`` within the function definitions.
* If ``normalize=true`` is set, the broadening factors are normalized to their values at :math:`L_\parallel=0` (i.e. upstream).

Output Variables
----------------

If the ``diagnose`` flag is set to ``true``, the component outputs the following variables:

* ``fieldline_geometry_parallel_length``: The parallel length along the flux tube.
* ``fieldline_geometry_geometric_broadening``: The geometric broadening factor.
* ``fieldline_geometry_transport_broadening``: The transport broadening factor.
* ``fieldline_geometry_flux_expansion``: The flux expansion factor.
* ``fieldline_geometry_flux_tube_broadening``: The total broadening factor applied to the flux tube.

Development status
----------------

This component is still being developed, and the interface is expected to change.