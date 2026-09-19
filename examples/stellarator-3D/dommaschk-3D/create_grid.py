import matplotlib.pyplot as plt
import numpy as np
import zoidberg
from zoidberg.fieldtracer import trace_poincare

# This is the matrix used for the Dommaschk potentials.
C = np.zeros((6, 5, 4))


# There are many different variations of Dommaschk potentials. We choose a quite simple version that already includes magnetic islands
C[5, 2, 1] = -1.489 / 1.0
C[5, 2, 2] = -1.489 / 1.0

# Magnetic field strength in Tesla
Btor = 1.5

# Major Radius in meter
R0 = 2.0

# This is the radial position of the inner flux surface at the outboard midplane
a1 = R0 + 0.08

# This is the radial position of the outer flux surface at the outboard midplane
a2 = R0 + 0.1

# n-fold symmetry. Because the coefficients are all 5-fold symmetric, we can simulate only one fifth of the whole device.
# symmetry = 1 results in a full-torus simulation
symmetry = 5

# Toroidal extend of the domain in rad
yperiod = 2.0 * np.pi / symmetry

# Grid size. Full torus results in ny = 40, while module simulations result in ny = 8
nx = 64 + 4
ny = int(5 * 8 / symmetry)
nz = 512

# Generate the actual magnetic field
field = zoidberg.field.DommaschkPotentials(C, R_0=R0, B_0=Btor)

# Toroidal coordinates
y_grid = np.linspace(0.0, yperiod, ny, endpoint=False)


from pathlib import Path

script_dir = Path(__file__).parent

fn = f"dommaschk_{nx}_{ny}_{nz}.fci.grid.nc"
filename = str(script_dir / fn)

plotting = False


# Don't create the grid when it is already existing
if (script_dir / fn).exists():
    print(
        filename,
        " exists. If you want to recreate the grid, please delete the currently existing one!",
    )

else:
    # This is the inner flux surface
    rzcoord, _ = trace_poincare(
        field,
        (a1),
        0.0,
        2.0 * np.pi / symmetry,
        y_slices=y_grid,
        revs=500,
        nplot=1,
        nover=20,
    )

    # This is the outer flux surface
    rzcoord2, _ = trace_poincare(
        field,
        (a2),
        0.0,
        2.0 * np.pi / symmetry,
        y_slices=y_grid,
        revs=500,
        nplot=1,
        nover=20,
    )

    # Here, the traced flux surfaces are converted into the lines required for Zoidberg
    inner_lines = []
    for i in range(ny):
        # spline_order = 1 is the most stable for most applications.
        inner_line = zoidberg.rzline.line_from_points(
            rzcoord[:, i, 0, 0], rzcoord[:, i, 0, 1], spline_order=1
        )
        # I find it more stable to downsample the lines to allways be lower than the actual grid size.
        # This would reduce grid-scale oscillations in the metric coefficients
        inner_line = inner_line.equallySpaced(n=nz // 4)
        inner_lines.append(inner_line)

    outer_lines = []
    for i in range(ny):
        # spline_order = 1 is the most stable for most applications.
        outer_line = zoidberg.rzline.line_from_points(
            rzcoord2[:, i, 0, 0], rzcoord2[:, i, 0, 1], spline_order=1
        )
        # I find it more stable to downsample the lines to allways be lower than the actual grid size.
        # This would reduce grid-scale oscillations in the metric coefficients
        outer_line = outer_line.equallySpaced(n=nz // 4)
        outer_lines.append(outer_line)

    if plotting:
        cs = 0
        fig, ax = plt.subplots()
        ax.scatter(
            rzcoord[:, cs, :, 0],
            rzcoord[:, cs, :, 1],
            edgecolor="None",
            color="black",
            s=2.0,
        )
        ax.set_aspect("equal")
        ax.grid(True)
        ax.set_ylim(-0.4, 0.4)
        ax.set_xlim(1.7, 2.3)
        plt.show()

    pol_grids = []

    for i in range(ny):
        # Using the elliptic grid generator to create the grid between the inner and outer flux surface
        pol_grid = zoidberg.poloidal_grid.grid_elliptic(
            inner_lines[i],
            outer_lines[i],
            nx,
            nz,
            restrict_size=2560,
            align=0,
            inner_ort=1,
            inner_maxmode=4,
            nx_inner=2,
            nx_outer=2,
        )
        pol_grids.append(pol_grid)

        if plotting:
            fig, ax = plt.subplots(figsize=(8, 8), dpi=400)
            pol_grid.plot(axis=ax, show=False)
            ax.set_aspect("equal")
            ax.grid(True, zorder=0)
            ax.set_xlim(1.6, 2.4)
            ax.set_ylim(-0.3, 0.3)
            plt.show()

    # y-coordinates of the grid slices
    ycoords = np.linspace(0, 2.0 * np.pi / symmetry, ny, endpoint=False)

    grid = zoidberg.grid.Grid(
        pol_grids, ycoords, 2.0 * np.pi / symmetry, yperiodic=True
    )

    # These are the maps of the tracing that contain the information for the parallel derivatives
    nslice = 1
    maps = zoidberg.make_maps(grid, field, nslice=nslice)

    # This step is crucial. The tracing is allways a little bit inaccurate. This results in some fieldlines exiting the flux surface
    # by incremental distances. This results in application of boundary conditions at all of these positions.
    # Here, we manually 'bend' back the first and last cells that are simulated back to allways be in the domain.
    maps["forward_xt_prime"][2, :, :] = 2.0
    maps["backward_xt_prime"][2, :, :] = 2.0
    maps["backward_xt_prime"][-3, :, :] = nx - 3
    maps["forward_xt_prime"][-3, :, :] = nx - 3

    if nslice == 2:
        maps["forward_xt_prime_2"][2, :, :] = 2.0
        maps["backward_xt_prime_2"][2, :, :] = 2.0
        maps["backward_xt_prime_2"][-3, :, :] = nx - 3
        maps["forward_xt_prime_2"][-3, :, :] = nx - 3

    # Writing the complete grid file
    with zoidberg.zoidberg.MapWriter(filename) as mw:
        mw.add_grid_field(grid, field)
        mw.add_maps(maps)
        mw.add_dagp()
    print("Finished creating the dommaschk grid!")
