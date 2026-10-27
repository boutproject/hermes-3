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
* ``fieldline_radius``: :math:`R`, the major radius of the magnetic fieldline
* ``poloidal_magnetic_field``: :math:`B_{pol}`, the poloidal magnetic field strength
* ``compute_Btor_from_R``: if ``true``, :math:`B_{tor} = B_{tor,up}R_{up}/R`, else :math:`B_{tor}` must be provided
* ``upstream_toroidal_magnetic_field``: :math:`B_{tor,up}`, required if ``compute_Btor_from_R=true``, the upstream toroidal magnetic field strength
* ``toroidal_magnetic_field``: :math:`B_{tor}`, required if ``compute_Btor_from_R=false``, the toroidal magnetic field along the fieldline
* ``diagnose``: add output variables to the state

For ``lambda_int``, ``fieldline_radius``, ``poloidal_magnetic_field`` and ``toroidal_magnetic_field``, you can provide expressions (as described in the `BOUT++ documentation <https://bout-dev.readthedocs.io/en/latest/user_docs/variable_init.html#expressions>`_).
In addition to the standard ``y`` and ``t`` expressions (``x`` and ``z`` don't vary in 1D), you can use ``{lpar}`` to write expressions in terms of the parallel length :math:`L_\parallel`.
For example,

.. code-block:: ini

 [fieldline_geometry]
 lambda_q = 1e-3
 b = 2.0
 lambda_int = where({lpar} > mesh:length_xpt, b, 1) * lambda_q

would set ``lambda_int`` to :math:`1mm` above the X-point and :math:`2mm` below the X-point.

Return variables
---------------------------

* ``fieldline_geometry_lpar``: :math:`L_\parallel`, parallel distance from upstream
* ``fieldline_geometry_lambda_int``: :math:`\lambda_{int}=\lambda_q + 1.64S`, the provided ``lambda_int`` function
* ``fieldline_geometry_magnetic_pitch``: :math:`\sin(\theta)=B_{pol}/B`, the magnetic pitch (ratio of poloidal to total field)
* ``fieldline_geometry_Rxy``: :math:`R`, the provided ``fieldline_radius`` function
* ``fieldline_geometry_Bpxy``: :math:`B_{pol}`, the poloidal magnetic field strength
* ``fieldline_geometry_Btxy``: :math:`B_{tor}`, the toroidal magnetic field strength
* ``fieldline_geometry_Bxy``: :math:`B`, the total magnetic field strength
* ``fieldline_geometry_transport_broadening``: :math:`\lambda_{int}/\lambda_{int,up}`, the flux tube broadening due to cross-field transport
* ``fieldline_geometry_f_R``: :math:`f_R=(B_{pol,up}/B_{up})/(B_{pol}/B)`, the total flux expansion
* ``fieldline_geometry_flux_tube_width``: :math:`\lambda_{int} f_R`, the flux tube radial width
* ``fieldline_geometry_dlpol``: :math:`dl_{pol}=dl_\parallel B_{pol}/B`, the poloidal length of the flux tube
* ``fieldline_geometry_cell_side_area``: :math:`2 \pi R \cdot dl_{pol}`, the poloidal length of the flux tube times its circumference
* ``fieldline_geometry_cell_volume``: :math:`2 \pi R \cdot dl_{pol} \lambda_{int} f_R`, the poloidal length of the flux tube times its circumference and radial width

Effect of these terms
---------------------------

The ``transport_broadening`` and ``magnetic_field_strength`` terms directly affect the Jacobian, which is set to :math:`J=(\lambda_{int}/\lambda_{q,up})/B`.
The ``cell_side_area`` and ``cell_volume`` affect the interaction of the cells with other components such as the reservoir model.


Reading inputs from a 1D grid file
----------------------------------

The magnetic geometry and transport broadening can also be read from a 1D grid file.
This is useful for a more precise comparison with a higher fidelity simulation, such as SOLPS-ITER, by using profiles extracted along a field line.

An example of a 1D grid file is shown below:

Example
^^^^^^^

This writes a grid file with the same mesh as ``examples/tokamak-1D/1D-threshold`` and a simple analytic geometry, however the process for using profiles taken from a higher fidelity simulation is the same:

.. code-block:: python

 import numpy as np
 from boututils.datafile import DataFile

 # Same mesh as examples/tokamak-1D/1D-threshold
 ny = 400
 length = 30.0  # [m]
 length_xpt = 10.0  # [m]
 dymin = 0.1

 y = 2 * np.pi * (np.arange(ny) + 0.5) / ny  # BOUT++ y at cell centres
 dy = (length / ny) * (1 + (1 - dymin) * (1 - y / np.pi))  # Cell widths [m]
 lpar = np.cumsum(dy) - 0.5 * dy  # Parallel distance of cell centres [m]

 # Simple geometry: R decreases and Bpol falls towards the target
 R = 1.4 - 0.4 * lpar / length  # [m]
 Bpol = 0.3 - 0.2 * lpar / length  # [T]
 Btor = -2.0 * R[0] / R  # [T]

 # 1 mm up to the X-point, broadening linearly to 3 mm at the target
 lambda_int = np.interp(lpar, [0, length_xpt, length], [1e-3, 1e-3, 3e-3])  # [m]

 with DataFile("grid_1d.nc", write=True, create=True) as f:
     f.write("nx", 1)
     f.write("ny", ny)
     f.write("nz", 1)
     f.write("ixseps1", -1)
     f.write("ixseps2", -1)
     f.write("dy", dy.reshape(1, ny))
     f.write("J", np.ones((1, ny)))
     f.write("lambda_int", lambda_int.reshape(1, ny))
     f.write("fieldline_radius", R.reshape(1, ny))
     f.write("poloidal_magnetic_field", Bpol.reshape(1, ny))
     f.write("toroidal_magnetic_field", Btor.reshape(1, ny))

Finally, point the ``[mesh]`` section at the grid file, and remove ``nx``, ``ny``, ``nz``, ``dy``, ``J``, ``ixseps1`` and ``ixseps2`` from ``[mesh]``, as these now come from the grid file:

.. code-block:: ini

 [mesh]
 file = "path/to/grid_1d.nc"

The ``fieldline_geometry`` component reads any of its inputs found in the grid file, which take precedence over the same options in the input file:

.. code-block:: ini

 [fieldline_geometry]
 # lambda_int, fieldline_radius, poloidal_magnetic_field and
 # toroidal_magnetic_field are all read from grid_1d.nc
 compute_Btor_from_R = false  # Use toroidal_magnetic_field from the grid file
 diagnose = true