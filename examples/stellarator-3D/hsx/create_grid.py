import fusionsc
import matplotlib.pyplot as plt
import numpy as np
import zoidberg

MXG = 2
ny = 10
nx = 256 + 2 * MXG
nz = 768


fscfield = fusionsc.devices.hsx.mainField(
    [-1.0722e04, -1.0722e04, -1.0722e04, -1.0722e04, -1.0722e04, -1.0722e04]
)
geometrygrid = fusionsc.devices.hsx.defaultGeometryGrid()
grid = fusionsc.devices.hsx.defaultGrid()
geometry = fusionsc.devices.hsx.wall()

fscfield = fscfield.compute(grid)
# %%
phiPlanes = np.linspace(0, 2 / 4 * np.pi, ny, endpoint=False)

# Put more starting points inside the island
rStart = np.linspace(1.2, 1.6, 30)
phiStart = 0.0

# Starting on an r, phi grid allows us to see all 5 islands
rStartG, phiStartG = np.meshgrid(rStart, phiStart, indexing="ij")

xStart = rStartG * np.cos(phiStartG)
yStart = rStartG * np.sin(phiStartG)
zStart = 0 * rStartG

startPoints = [xStart, yStart, zStart]


# %%

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


# %%


inner_lines = []
outer_lines = []
pol_grids = []


# iPhi = 8

for iPhi in range(ny):
    fs_inner = 15
    fs_axis = 18

    x_inner, y_inner, z_inner, _, _ = pcPoints[:, iPhi, fs_inner, :]
    r_inner = np.sqrt(x_inner**2 + y_inner**2)

    x_axis, y_axis, z_axis, _, _ = pcPoints[:, iPhi, fs_axis, :]
    r_axis = np.sqrt(x_axis**2 + y_axis**2)

    r_center = np.mean(r_axis)
    z_center = np.mean(z_axis)

    x_wall, y_wall, z_wall = geometry.planarCut(phi=phiPlanes[iPhi])
    r_wall = np.sqrt(x_wall**2 + y_wall**2)

    # fig,ax = plt.subplots()
    # ax.scatter(r_inner[0], z_inner[0])

    # ax.scatter(r_wall[::10,0], z_wall[::10,0])

    # plt.show()

    inner = zoidberg.rzline.line_from_points(
        r_inner[0], z_inner[0], spline_order=1
    ).equallySpaced(int(2 * nz))
    outer = zoidberg.rzline.line_from_points(
        r_wall[::5, 0], z_wall[::5, 0], spline_order=1
    ).equallySpaced(int(2 * nz))

    inner_lines.append(inner)
    outer_lines.append(outer)

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

    fig, ax = plt.subplots()
    # inner.plot(axis=ax,show = False)
    # outer.plot(axis=ax,show = False)
    ax.set_aspect("equal")
    pol_grid.plot(axis=ax)


# %%
symmetry = 4
ycoords = np.linspace(0, 2.0 * np.pi / symmetry, ny, endpoint=False)
magnetic_field = zoidberg.field.FusionSCField(fscfield)
zbgrid = zoidberg.grid.Grid(pol_grids, ycoords, 2.0 * np.pi / symmetry, yperiodic=True)
maps = zoidberg.make_maps(zbgrid, magnetic_field, nslice=1)
maps["forward_xt_prime"][2, :, :] = 2.0
maps["backward_xt_prime"][2, :, :] = 2.0
fn = f"hsx_{nx}_{ny}_{nz}.fci.grid.nc"

zoidberg.zoidberg.write_maps(zbgrid, magnetic_field, maps, metric2d=False, gridfile=fn)
