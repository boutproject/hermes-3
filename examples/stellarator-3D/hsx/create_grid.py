import fusionsc
import numpy as np
import zoidberg

# Number of guard cells
MXG = 2

# Number of cross sections
ny = 10

# Radial cells
nx = 256 + 2 * MXG

# Poloidal cells
nz = 768

# Coil currents
fscfield = fusionsc.devices.hsx.mainField(
    [-1.0722e04, -1.0722e04, -1.0722e04, -1.0722e04, -1.0722e04, -1.0722e04]
)
# Use all default FusionSC components
geometrygrid = fusionsc.devices.hsx.defaultGeometryGrid()
grid = fusionsc.devices.hsx.defaultGrid()
geometry = fusionsc.devices.hsx.wall()
fscfield = fscfield.compute(grid)

# Symmetry: 4 = one module, 1 = full torus
symmetry = 4
phiPlanes = np.linspace(0, 2 * np.pi / symmetry, ny, endpoint=False)

# Put more starting points inside the island
rStart = np.linspace(1.2, 1.6, 30)
phiStart = 0.0

# Starting on an r, phi grid allows us to see all 5 islands
rStartG, phiStartG = np.meshgrid(rStart, phiStart, indexing="ij")

xStart = rStartG * np.cos(phiStartG)
yStart = rStartG * np.sin(phiStartG)
zStart = 0 * rStartG

startPoints = [xStart, yStart, zStart]

# Tracing by FusionSC
pcPoints = fusionsc.flt.poincareInPhiPlanes(
    startPoints,
    fscfield,
    phiPlanes,
    512,
    grid=grid,
    geometry=geometry,
    geometryGrid=geometrygrid,
    distanceLimit=1e4,
    targetError=1e-5,
    maxStepSize=0.1,
)


inner_lines = []
outer_lines = []
pol_grids = []

# For each cross section
for iPhi in range(ny):
    # This is the flux surface of the poincare plots that will
    # be the inner contour of the grid
    fs_inner = 15

    # This is a flux surface close to the axis
    fs_axis = 18

    x_inner, y_inner, z_inner, _, _ = pcPoints[:, iPhi, fs_inner, :]
    r_inner = np.sqrt(x_inner**2 + y_inner**2)

    x_axis, y_axis, z_axis, _, _ = pcPoints[:, iPhi, fs_axis, :]
    r_axis = np.sqrt(x_axis**2 + y_axis**2)

    r_center = np.mean(r_axis)
    z_center = np.mean(z_axis)

    x_wall, y_wall, z_wall = geometry.planarCut(phi=phiPlanes[iPhi])
    r_wall = np.sqrt(x_wall**2 + y_wall**2)

    # Create the inner and outer contour from the poincare plot
    # spline_order = 1 helps with bad polynomial fits
    # equallySpaced(2*nz) is a choice I usually do
    # The grid generator seems like prefer when the inner/outer lines
    # Have more points than the nz of the grid
    inner = zoidberg.rzline.line_from_points(
        r_inner[0], z_inner[0], spline_order=1
    ).equallySpaced(int(2 * nz))
    outer = zoidberg.rzline.line_from_points(
        r_wall[::5, 0], z_wall[::5, 0], spline_order=1
    ).equallySpaced(int(2 * nz))

    inner_lines.append(inner)
    outer_lines.append(outer)

    # Elliptic grid generator
    pol_grid = zoidberg.poloidal_grid.grid_elliptic(
        inner,
        outer,
        nx,
        nz // 2,
        inner_maxmode=1,
        dz_relax=1.75,
        nx_inner=2,
        nx_outer=2,
    )

    pol_grids.append(pol_grid)


# Toroidal coordinates of each cross section
ycoords = np.linspace(0, 2.0 * np.pi / symmetry, ny, endpoint=False)

# Magnetic field, now in Zoidberg
magnetic_field = zoidberg.field.FusionSCField(fscfield)

# The grid that Zoidberg needs for the full gridfile
zbgrid = zoidberg.grid.Grid(pol_grids, ycoords, 2.0 * np.pi / symmetry, yperiodic=True)

# This does the tracing, usually the most expensive part
maps = zoidberg.make_maps(zbgrid, magnetic_field, nslice=1)

# This adjust the most inner ring of cells to land exactly on themselves in index-space
# It prevents possible application of boundary conditions, when field lines
# enter the core. This should in principle not happen, but finite tolerances lead
# to shallow boundary crossings

maps["forward_xt_prime"][2, :, :] = 2.0
maps["backward_xt_prime"][2, :, :] = 2.0

fn = f"hsx_{nx}_{ny}_{nz}.fci.grid.nc"

# Final step, writing the grid
zoidberg.zoidberg.write_maps(zbgrid, magnetic_field, maps, metric2d=False, gridfile=fn)
