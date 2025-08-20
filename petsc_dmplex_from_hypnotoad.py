from plot_corners_functions import plot_corners_get_dmplex_data
from petsc4py import PETSc
import meshio

# function to check correctly assigned labels
def check_label_value_coords(value,coords,boundary_vertex_info):
    check = False
    for key in boundary_vertex_info.keys():
        if value == boundary_vertex_info[key]["DMFaceSetsLabel"]:
            Rxy = boundary_vertex_info[key]["Rxy"]
            Zxy = boundary_vertex_info[key]["Zxy"]
            if coords[0] in Rxy and coords[1] in Zxy:
                check = True
                break
    return check
def create_and_visualize_mesh(file_path,Nx,Ny,cell_list,vertex_list,boundary_vertex_info):
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
    label_name = "Face Sets"
    dm.createLabel(label_name)
    global_label = False
    if global_label:
        dm.markBoundaryFaces(label_name, value=100)
    else:
        # physical labels
        # Get depth integer to permit extracting all faces (depth - 1)
        depth = dm.getDepth()
        # get min/max indices for faces
        face_start, face_end = dm.getDepthStratum(depth - 1)
        # get min/max indices for vertices
        vertex_start, vertex_end = dm.getDepthStratum(depth - 2)
        # loop over all boundary types
        for key in boundary_vertex_info.keys():
            label_value = boundary_vertex_info[key]["DMFaceSetsLabel"]
            boundary_vertex_indices = boundary_vertex_info[key]["ivertex"]
            print(key)
            print(boundary_vertex_indices)
            # loop over faces
            for face in range(face_start, face_end):
                # vertices supporting this face
                cone = dm.getCone(face)
                # n.b. need to subtract vertex_start
                # to get back the input global index
                # if all points in cone are in the boundary, this is a boundary face
                if all( (v-vertex_start) in boundary_vertex_indices for v in cone):
                    dm.setLabelValue("Face Sets", face, label_value)

    # Print labeled boundary faces
    depth = dm.getDepth()
    face_start, face_end = dm.getDepthStratum(depth - 1)
    vertex_start, vertex_end = dm.getDepthStratum(depth - 2)
    print(f"Boundary faces labeled with {label_name}:")
    for face in range(face_start, face_end):
        value = dm.getLabelValue(label_name, face)
        cone = dm.getCone(face)
        coords = vertex_coords[cone - vertex_start]
        if value in [100,200,300,400]:
            print(f"Face {face} labeled with value {value}")
            print(f"Face {face} coords {coords[0]}")
            check = check_label_value_coords(value,coords[0],boundary_vertex_info)
            print(f"Test passed: {check}")
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
Nx,Ny,cell_vertices,vertex_list,boundary_vertex_info = plot_corners_get_dmplex_data(file_path)
create_and_visualize_mesh(file_path,Nx,Ny,cell_vertices,vertex_list,boundary_vertex_info)
