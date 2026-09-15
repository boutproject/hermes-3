import matplotlib.pyplot as plt
import numpy as np
import zoidberg
from zoidberg.fieldtracer import trace_poincare

C = np.zeros((6, 5, 4))

C[5, 2, 1] = -1.489 / 1.0
C[5, 2, 2] = -1.489 / 1.0

Btor = 1.5
R0 = 2.0
a1 = R0 + 0.08
a2 = R0 + 0.1

symmetry = 5.0
yperiod = 2.0 * np.pi / symmetry

nx = 32 + 4
ny = int(5 * 8 / symmetry)
nz = 256

field = zoidberg.field.DommaschkPotentials(C, R_0=R0, B_0=Btor)

y_grid = np.linspace(0.0, yperiod, ny, endpoint=False)


from pathlib import Path

script_dir = Path(__file__).parent

fn = f"dommaschk_{nx}_{ny}_{nz}.fci.grid.nc"
filename = str(script_dir / fn)

plotting = False

if (script_dir / fn).exists():
    print(
        filename,
        " exists. If you want to recreate the grid, please delete the currently existing one!",
    )

else:
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

    inner_lines = []
    for i in range(ny):
        inner_line = zoidberg.rzline.line_from_points(
            rzcoord[:, i, 0, 0], rzcoord[:, i, 0, 1], spline_order=1
        )
        inner_line = inner_line.equallySpaced(n=nz // 4)
        inner_lines.append(inner_line)

    outer_lines = []
    for i in range(ny):
        outer_line = zoidberg.rzline.line_from_points(
            rzcoord2[:, i, 0, 0], rzcoord2[:, i, 0, 1], spline_order=1
        )
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

    ycoords = np.linspace(0, 2.0 * np.pi / symmetry, ny, endpoint=False)

    grid = zoidberg.grid.Grid(
        pol_grids, ycoords, 2.0 * np.pi / symmetry, yperiodic=True
    )

    maps = zoidberg.make_maps(grid, field, nslice=2)

    maps["forward_xt_prime"][2, :, :] = 2.0
    maps["backward_xt_prime"][2, :, :] = 2.0
    maps["backward_xt_prime"][-3, :, :] = nx - 3
    maps["forward_xt_prime"][-3, :, :] = nx - 3

    maps["forward_xt_prime_2"][2, :, :] = 2.0
    maps["backward_xt_prime_2"][2, :, :] = 2.0
    maps["backward_xt_prime_2"][-3, :, :] = nx - 3
    maps["forward_xt_prime_2"][-3, :, :] = nx - 3

    zoidberg.write_maps(grid, field, maps, metric2d=False, gridfile=filename)

    print("Finished creating the dommaschk grid!")
