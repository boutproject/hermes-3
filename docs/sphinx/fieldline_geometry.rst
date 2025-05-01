.. _sec-fieldline_geometry:

Fieldline Geometry Component
============================

The ``fieldline_geometry`` component in Hermes-3 helps to set up 1D simulations and to compute cell geometry for other modules.

Implementation
---------------------------

The component is implemented in ``hermes-3/src/fieldline_geometry.cxx`` and its header file is ``hermes-3/include/fieldline_geometry.hxx``.

Input parameters
---------------------------

* ``lambda_int``: :math:`\lambda_q + 1.64S`, the radial width of the flux tube mapped upstream (i.e. removing flux expansion effects)
* ``fieldline_pitch_angle``: :math:`\sin(\theta)=B_{pol}/B`, the pitch angle of the magnetic field
* ``fieldline_radius``: :math:`R`, the major radius of the magnetic fieldline
* ``compute_B_from_R``: if ``true``, :math:`B=B_{up}R_{up}/R`, else :math:`B` must be provided
* ``upstream_magnetic_field_strength``: :math:`B_{up}`, required if ``compute_B_from_R=true``, the upstream magnetic field strength
* ``magnetic_field_strength``: :math:`B`, required if ``compute_B_from_R=false``, the magnetic field along the fieldline
* ``diagnose``: add output variables to the state

For ``lambda_int``, ``fieldline_pitch_angle``, ``fieldline_radius`` and ``magnetic_field_strength``, you can provide expressions (as described in the `BOUT++ documentation <https://bout-dev.readthedocs.io/en/latest/user_docs/variable_init.html#expressions>`_). In addition to the standard ``y`` and ``t`` expressions (``x`` and ``z`` don't vary in 1D), you can use ``{lpar}`` to write expressions in terms of the parallel length :math:`L_\parallel`. For example,

.. code-block:: ini

 [fieldline_geometry]
 lambda_q = 1e-3
 b = 2.0
 lambda_int = where({lpar} > mesh:length_xpt, 1, b) * lambda_q

would set ``lambda_int`` to :math:`1mm` above the X-point and :math:`2mm` below the X-point.

Return variables
---------------------------

* ``fieldline_geometry_lpar``: :math:`L_\parallel`, parallel distance from upstream
* ``fieldline_geometry_lambda_int``: :math:`\lambda_{int}=\lambda_q + 1.64S`, the provided ``lambda_int`` function
* ``fieldline_geometry_pitch_angle``: :math:`\sin(\theta)=B_{pol}/B`, the provided ``fieldline_pitch_angle`` function
* ``fieldline_geometry_fieldline_radius``: :math:`R`, the provided ``fieldline_radius`` function
* ``fieldline_geometry_magnetic_field_strength``: :math:`B`, the provided or computed ``magnetic_field_strength`` function
* ``fieldline_geometry_transport_broadening``: :math:`\lambda_{int}/\lambda_{q,up}`, the flux tube broadening due to cross-field transport
* ``fieldline_geometry_flux_expansion``: :math:`f_{exp}=(B_{pol,up}/B_{up})/(B_{pol}/B)`, the flux tube broadening due to flux expansion
* ``fieldline_geometry_flux_tube_width``: :math:`\lambda_{int} f_{exp}`, the flux tube radial width
* ``fieldline_geometry_cell_poloidal_length``: :math:`dl_{pol}=dl_\parallel B_{pol}/B`, the poloidal length of the flux tube
* ``fieldline_geometry_cell_side_area``: :math:`2 \pi R \cdot dl_{pol}`, the poloidal length of the flux tube times its circumference
* ``fieldline_geometry_cell_volume``: :math:`2 \pi R \cdot dl_{pol} \lambda_{int} f_{exp}`, the poloidal length of the flux tube times its circumference and radial width

Effect of these terms
---------------------------

The ``transport_broadening`` and ``magnetic_field_strength`` terms directly affect the Jacobian, which is set to :math:`J=(\lambda_{int}/\lambda_{q,up})/B`.

The ``cell_side_area`` and ``cell_volume`` affect the interaction of the cells with other components such as the reservoir model.

