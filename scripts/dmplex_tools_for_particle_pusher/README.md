DMPlex tools for pushing particles
==================================

The tools in this directory are designed to work with the [DMPlex](https://petsc.org/release/manualpages/DMPlex/)
mesh format. This format is used by [NESO-Particles](https://github.com/ExCALIBUR-NEPTUNE/NESO-Particles)
to define the mesh interface for codes using NESO-Particles as a library.

In Hermes-3, with [BOUT++](https://github.com/boutproject/BOUT-dev) we use [Hypnotoad](https://github.com/boutproject/hypnotoad) to generate a NetCDF file with grid data for the conservative finite-difference
scheme used by Hermes-3.

In exploratory work we created `plot_corners_functions.py` and `find_global_vertices_save_to_hypnotoad_nc.py` which can be used to get the necessary data from the Hypnotoad grid file to make a global list of corners (vertices) of the Hypnotoad grid, and assemble the data in the right format to give to PETSc to create a DMPlex. The scripts can also save a
modified version of the Hypnotoad grid to include this data. Note that some of the functionality here is now captured by `../../bout-particle-push.cxx`.

The following scripts use `argparse` to take command line arguments.
* `find_global_vertices_save_to_hypnotoad_nc.py`: find the global vertices list of a Hynotoad grid file and save the results to a modified Hypnotoad file. Plot the global vertices identified.
* `petsc_dmplex_from_hypnotoad.py`: this python script uses [petsc4py](https://petsc.org/main/petsc4py/reference/petsc4py.PETSc.html) to create a `.h5` file that represents the DMPlex data for a mesh that approximates a Hypnotoad grid.
* `load_and_plot_dm.py`: this program script the edges defining the DMPlex grid from a `.h5` file created by PETSc.
* `particle_animator`: this program script the DMPlex grid using the same methods as `load_and_plot_dm.py` and then plots particle trajectories created using the `../../bout-particle-push.cxx` main program.