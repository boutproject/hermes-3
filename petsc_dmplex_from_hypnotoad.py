from plot_corners_functions import plot_corners_get_dmplex_data
from petsc4py import PETSc
import meshio

def create_and_visualize_mesh(file_path,Nx,Ny,cell_list,vertex_list):
    # Initialize PETSc
    comm = PETSc.COMM_WORLD
    # Create a DMPlex object
    dm = PETSc.DMPlex().create(comm=comm)
    # Define mesh parameters
    # for a rectangle with 2 cells and 6 vertices
    dim = 2  # 2D mesh
    num_cells = Nx*Ny
    # integers defining the cells
    #cells = [[0, 1, 3, 4], [1, 2, 4, 5]]
    cells = cell_list.astype('int32')
    print("cells")
    print(cells)
    # global list of vertices, in correct order
    # vertex_coords = [[0.0, 0.0],
    #                  [0.0, 1.0],
    #                  [1.0, 0.0],
    #                  [1.0, 1.0],
    #                  [2.0, 0.0],
    #                  [2.0, 1.0],
    #                 ]
    # this is the incorrect order to supply the vertices
    # vertex_coords = [[0.0, 0.0],
    #                  [1.0, 0.0],
    #                  [2.0, 0.0],
    #                  [0.0, 1.0],
    #                  [1.0, 1.0],
    #                  [2.0, 1.0],
    #                 ]
    vertex_coords = vertex_list
    print("vertex_coords")
    print(vertex_coords)
    # Create the mesh from the cell list
    dm.createFromCellList(dim, cells, vertex_coords, comm=comm)
    # Distribute the mesh (optional, useful for parallel runs)
    #dm = dm.distribute()
    mesh_name = file_path+".mesh"
    # Write the mesh to a VTK file
    viewer = PETSc.Viewer().createVTK(mesh_name+".vtk", mode=PETSc.Viewer.Mode.WRITE, comm=comm)
    dm.view(viewer)
    viewer.destroy()
    # Write the mesh to a .h5 file
    #viewer = PETSc.Viewer().createHDF5('mesh.h5', mode=PETSc.Viewer.Mode.WRITE, comm=comm)
    viewer = PETSc.ViewerHDF5().create(mesh_name+'.h5', mode=PETSc.Viewer.Mode.WRITE, comm=comm)
    dm.view(viewer)
    viewer.destroy()
    # Read the mesh from .h5
    dmtest = PETSc.DMPlex().create(comm=comm)
    viewer = PETSc.ViewerHDF5().create(mesh_name+'.h5', mode=PETSc.Viewer.Mode.READ, comm=comm)
    dmtest.load(viewer)
    viewer.destroy()
    # Write the reloaded mesh to a VTK file
    viewer = PETSc.Viewer().createVTK(mesh_name+"_h5_to_vtk.vtk", mode=PETSc.Viewer.Mode.WRITE, comm=comm)
    dmtest.view(viewer)
    viewer.destroy()
    # write mesh to gmsh format with ascii output, for use in NESO-particles tests
    mesh = meshio.read(mesh_name+".vtk")
    mesh.write(mesh_name+".msh",binary=False)
    return None

file_path = 'dmtest_data/expected_nonorthogonal.grd.nc'
#file_path = 'dmtest_data/nonorthog.bout.grd.nc'
#file_path = 'dmtest_data/expected_orthogonal.grd.nc'
#file_path = 'dmtest_data/udn.bout.grd.nc'
#file_path = 'dmtest_data/lsn.bout.grd.nc'
Nx,Ny,cell_vertices,vertex_list = plot_corners_get_dmplex_data(file_path)
create_and_visualize_mesh(file_path,Nx,Ny,cell_vertices,vertex_list)
