import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import os

def save_ivertex_indices_to_netcdf(source_file,destination_file,Rpoints_full,Zpoints_full):
    # copy source file to destination
    os.system(f"cp {source_file} {destination_file}")
    # open destination file
    dataset = nc.Dataset(destination_file, mode="a")

    # get data including y guards
    # but note that Rpoints_full,Zpoints_full created
    # without y guards as the DMPlex should not have guards
    data_Rxy = np.copy(dataset.variables['Rxy_corners'][:])
    data_Zxy = np.copy(dataset.variables['Zxy_corners'][:])
    data_Rxy_lr = np.copy(dataset.variables['Rxy_lower_right_corners'][:])
    data_Zxy_lr = np.copy(dataset.variables['Zxy_lower_right_corners'][:])
    data_Rxy_ur = np.copy(dataset.variables['Rxy_upper_right_corners'][:])
    data_Zxy_ur = np.copy(dataset.variables['Zxy_upper_right_corners'][:])
    data_Rxy_ul = np.copy(dataset.variables['Rxy_upper_left_corners'][:])
    data_Zxy_ul = np.copy(dataset.variables['Zxy_upper_left_corners'][:])

    ivertex_corners = get_boutxx_corner_index(data_Rxy,data_Zxy,Rpoints_full,Zpoints_full)
    ivertex_corners_lr = get_boutxx_corner_index(data_Rxy_lr,data_Zxy_lr,Rpoints_full,Zpoints_full)
    ivertex_corners_ul = get_boutxx_corner_index(data_Rxy_ul,data_Zxy_ul,Rpoints_full,Zpoints_full)
    ivertex_corners_ur = get_boutxx_corner_index(data_Rxy_ur,data_Zxy_ur,Rpoints_full,Zpoints_full)
    # index saved as double as BOUT++ can only handle double Field2D
    ptr = dataset.createVariable("ivertex_lower_left_corners","f8",("x", "y"))
    ptr.setncattr('bout_type','Field2D')
    ptr[:] = ivertex_corners
    ptr = dataset.createVariable("ivertex_lower_right_corners","f8",("x", "y"))
    ptr.setncattr('bout_type','Field2D')
    ptr[:] = ivertex_corners_lr
    ptr = dataset.createVariable("ivertex_upper_right_corners","f8",("x", "y"))
    ptr.setncattr('bout_type','Field2D')
    ptr[:] = ivertex_corners_ur
    ptr = dataset.createVariable("ivertex_upper_left_corners","f8",("x", "y"))
    ptr.setncattr('bout_type','Field2D')
    ptr[:] = ivertex_corners_ul
    dataset.close()
    return None

def isapprox(a,b,tol=1.0e-8):
        if abs(a - b) < tol:
            return True
        else:
            return False

def unique_points_2D(corners_Rxy,corners_Zxy):
    Nx, Ny = corners_Rxy.shape
    unique=True
    for ix in range(0,Nx):
        for iy in range(0,Ny):
            for ixp in range(0,Nx):
                for iyp in range(0,Ny):
                    if ix == ixp and iy == iyp:
                        continue
                    if (isapprox(corners_Rxy[ix,iy],corners_Rxy[ixp,iyp]) and isapprox(corners_Zxy[ix,iy],corners_Zxy[ixp,iyp])):
                        ic = iy + Ny*ix
                        print("non-unique point:",ix,iy,ic)
                        unique = False
    print("unique points =",unique)
    return None

def unique_points_1D(corners_Rxy,corners_Zxy):
    Nc, = corners_Rxy.shape
    unique=True
    for ic in range(0,Nc):
        for icp in range(0,Nc):
            if ic == icp:
                continue
            if (isapprox(corners_Rxy[ic],corners_Rxy[icp]) and isapprox(corners_Zxy[ic],corners_Zxy[icp])):
                print("non-unique point:",ix,iy,ic)
                unique = False
    print("unique points =",unique)
    return None

def remove_nonunique_points(corners_Rxy,corners_Zxy):
    Nc, = corners_Rxy.shape
    duplicate_indices = []
    for ic in range(0,Nc):
        for icp in range(0,Nc):
            if ic == icp:
                continue
            if (isapprox(corners_Rxy[ic],corners_Rxy[icp]) and isapprox(corners_Zxy[ic],corners_Zxy[icp])):
                print("non-unique point:",icp)
                # store the duplicate index, if we have not already stored the index for this point
                if (ic not in duplicate_indices) and (icp not in duplicate_indices):
                    duplicate_indices.append(icp)
    print(duplicate_indices)
    if len(duplicate_indices) > 0:
        corners_Rxy = np.delete(corners_Rxy,duplicate_indices)
        corners_Zxy = np.delete(corners_Zxy,duplicate_indices)
    return corners_Rxy, corners_Zxy

# identify which vertex an (R,Z) pair, or print error message
def convert_R_Z_to_vertex_index(R,Z,R_vertices,Z_vertices,failure=None):
    Nc = len(R_vertices)
    for ic in range(0,Nc):
        if isapprox(R,R_vertices[ic]) and isapprox(Z,Z_vertices[ic]):
            ic_vertex = ic
            return ic_vertex
    print("Failed to find (R,Z) index")
    return failure

def get_cell_vertex_list_inner(R_ll,Z_ll, R_lr, Z_lr,
                        R_ur, Z_ur, R_ul, Z_ul,
                        R_vertices, Z_vertices):
    # R coords defining the vertcies
    R_local_vertex_list = np.array([R_ll, R_lr, R_ur, R_ul])
    # R location of cell centre
    R_mid = np.mean(R_local_vertex_list)
    # Z coords defining the vertices, in same order as R above
    Z_local_vertex_list = np.array([Z_ll, Z_lr, Z_ur, Z_ul])
    # Z location of cell centre
    Z_mid = np.mean(Z_local_vertex_list)
    # theta is the angle from the vertex to the
    # midpoint of the cell, measured
    # with respect to the major radial direction
    theta_list = np.zeros(4)
    for j in range(0,4):
        R = R_local_vertex_list[j] - R_mid
        Z = Z_local_vertex_list[j] - Z_mid
        theta_list[j] = np.arctan2(Z,R)
    # order theta in increasing order
    # to get the order for anticlockwise
    # specification of the vertices
    sort_indices = np.argsort(theta_list)
    R_local_vertex_sorted = R_local_vertex_list[sort_indices]
    Z_local_vertex_sorted = Z_local_vertex_list[sort_indices]
    vertices = []
    # append the vertices in anti-clockwise order
    for j in range(0,4):
        R = R_local_vertex_sorted[j]
        Z = Z_local_vertex_sorted[j]
        vertices.append(convert_R_Z_to_vertex_index(R,Z,R_vertices,Z_vertices))
    return vertices

def get_cell_vertex_list(Rxy_ll, Zxy_ll,
        Rxy_lr, Zxy_lr, Rxy_ur, Zxy_ur, Rxy_ul, Zxy_ul,
        R_vertices, Z_vertices):
    Nx, Ny = np.shape(Rxy_ll)
    Ncells = Nx*Ny
    cell_vertices = np.zeros((Ncells,4),dtype=int)
    cell_vertices_RZ = np.zeros((Ncells,4,2),dtype=float)
    cell_vertices[:,:] = -3000
    for ix in range(0,Nx):
        for iy in range(0,Ny):
            icell = iy + Ny*ix
            R_ll = Rxy_ll[ix,iy]
            Z_ll = Zxy_ll[ix,iy]
            R_lr = Rxy_lr[ix,iy]
            Z_lr = Zxy_lr[ix,iy]
            R_ur = Rxy_ur[ix,iy]
            Z_ur = Zxy_ur[ix,iy]
            R_ul = Rxy_ul[ix,iy]
            Z_ul = Zxy_ul[ix,iy]
            cell_vertices[icell,:] = get_cell_vertex_list_inner(R_ll,Z_ll, R_lr, Z_lr,
                        R_ur, Z_ur, R_ul, Z_ul,
                        R_vertices, Z_vertices)
            # this list must be in the same order as the list of cell vertices above
            cell_vertices_RZ[icell,0,:] = [R_ll,Z_ll]
            cell_vertices_RZ[icell,1,:] = [R_lr,Z_lr]
            cell_vertices_RZ[icell,2,:] = [R_ur,Z_ur]
            cell_vertices_RZ[icell,3,:] = [R_ul,Z_ul]
    # get a global list of vertices
    Nc, = np.shape(R_vertices)
    vertex_list = np.zeros((Nc,2),dtype=float)
    for ic in range(0,Nc):
        vertex_list[ic,0] = R_vertices[ic]
        vertex_list[ic,1] = Z_vertices[ic]
    return cell_vertices, cell_vertices_RZ, vertex_list

def get_boutxx_corner_index(Rxy_corners,Zxy_corners,R_vertices,Z_vertices):
    Nx, Ny = Rxy_corners.shape
    iglobal_corners = np.zeros((Nx,Ny),dtype=int)
    if Nx*Ny > len(R_vertices):
        failure_return = -1
    else:
        failure_return = None
    for ix in range(0,Nx):
        for iy in range(0,Ny):
            # find the global index for this corner, and store
            iglobal_corners[ix,iy] = convert_R_Z_to_vertex_index(Rxy_corners[ix,iy],Zxy_corners[ix,iy],
                                                                R_vertices,Z_vertices,failure=failure_return)
    return iglobal_corners

def get_pfr_lower_boundary_vertices(ivertex_corners,data_Rxy,data_Zxy,
                    ivertex_corners_ul,data_Rxy_ul,data_Zxy_ul,
                    y_boundary_guards,jyseps1_1,exclude_y_guard_cells=False):
    Nx, Ny = ivertex_corners.shape
    ivertex_pfr_lower = []
    Rxy_pfr_lower = []
    Zxy_pfr_lower = []
    if exclude_y_guard_cells:
        jy_bndry = 0
    else:
        jy_bndry = y_boundary_guards
    lim1 = jyseps1_1+jy_bndry+1
    for j in range(jy_bndry,lim1):
        ivertex_pfr_lower.append(ivertex_corners[0,j])
        Rxy_pfr_lower.append(data_Rxy[0,j])
        Zxy_pfr_lower.append(data_Zxy[0,j])
    lim2 = Ny - jy_bndry - jyseps1_1 - 1
    lim3 = Ny - jy_bndry
    for j in range(lim2,lim3):
        ivertex_pfr_lower.append(ivertex_corners[0,j])
        Rxy_pfr_lower.append(data_Rxy[0,j])
        Zxy_pfr_lower.append(data_Zxy[0,j])
    ivertex_pfr_lower.append(ivertex_corners_ul[0,lim3-1])
    Rxy_pfr_lower.append(data_Rxy_ul[0,lim3-1])
    Zxy_pfr_lower.append(data_Zxy_ul[0,lim3-1])
    return ivertex_pfr_lower, Rxy_pfr_lower, Zxy_pfr_lower

def get_pfr_upper_boundary_vertices(ivertex_corners,data_Rxy,data_Zxy,
                        ivertex_corners_ul,data_Rxy_ul,data_Zxy_ul,
                        y_boundary_guards,jyseps1_2,jyseps2_1,jyseps2_2,
                        ny_inner, exclude_y_guard_cells=False):
    ivertex_pfr_upper = []
    Rxy_pfr_upper = []
    Zxy_pfr_upper = []
    if exclude_y_guard_cells:
        jyseps2_1g = jyseps2_1
    else:
        jyseps2_1g = jyseps2_1 + y_boundary_guards
    for j in range(ny_inner,jyseps1_2+1):
        ivertex_pfr_upper.append(ivertex_corners[0,j])
        Rxy_pfr_upper.append(data_Rxy[0,j])
        Zxy_pfr_upper.append(data_Zxy[0,j])
    for j in range(jyseps2_1+1,ny_inner):
        ivertex_pfr_upper.append(ivertex_corners[0,j])
        Rxy_pfr_upper.append(data_Rxy[0,j])
        Zxy_pfr_upper.append(data_Zxy[0,j])
    j = ny_inner-1
    ivertex_pfr_upper.append(ivertex_corners_ul[0,j])
    Rxy_pfr_upper.append(data_Rxy_ul[0,j])
    Zxy_pfr_upper.append(data_Zxy_ul[0,j])
    return ivertex_pfr_upper, Rxy_pfr_upper, Zxy_pfr_upper

def get_vac_left_boundary_vertices(ivertex_corners_lr,data_Rxy_lr,data_Zxy_lr,
                                    ivertex_corners_ur,data_Rxy_ur,data_Zxy_ur,
                                    y_boundary_guards,ny_inner):
    ivertex_sol_vac_left = []
    Rxy_sol_vac_left = []
    Zxy_sol_vac_left = []
    jlim = ny_inner#+2*y_boundary_guards
    for j in range(0,jlim):
        ivertex_sol_vac_left.append(ivertex_corners_lr[-1,j])
        Rxy_sol_vac_left.append(data_Rxy_lr[-1,j])
        Zxy_sol_vac_left.append(data_Zxy_lr[-1,j])
    j = jlim - 1
    ivertex_sol_vac_left.append(ivertex_corners_ur[-1,j])
    Rxy_sol_vac_left.append(data_Rxy_ur[-1,j])
    Zxy_sol_vac_left.append(data_Zxy_ur[-1,j])
    return ivertex_sol_vac_left, Rxy_sol_vac_left, Zxy_sol_vac_left

def get_vac_right_boundary_vertices(ivertex_corners_lr,data_Rxy_lr,data_Zxy_lr,
                                    ivertex_corners_ur,data_Rxy_ur,data_Zxy_ur,
                                    jyseps1_2,ny_inner):
    Nx, Ny = ivertex_corners_lr.shape
    ivertex_sol_vac_right = []
    Rxy_sol_vac_right = []
    Zxy_sol_vac_right = []
    j = ny_inner#jyseps1_2+1
    ivertex_sol_vac_right.append(ivertex_corners_lr[-1,j])
    Rxy_sol_vac_right.append(data_Rxy_lr[-1,j])
    Zxy_sol_vac_right.append(data_Zxy_lr[-1,j])
    for j in range(ny_inner,Ny):#jyseps1_2+1
        ivertex_sol_vac_right.append(ivertex_corners_ur[-1,j])
        Rxy_sol_vac_right.append(data_Rxy_ur[-1,j])
        Zxy_sol_vac_right.append(data_Zxy_ur[-1,j])
    return ivertex_sol_vac_right, Rxy_sol_vac_right, Zxy_sol_vac_right

def get_core_boundary_vertices(ivertex_corners, data_Rxy, data_Zxy,
   jyseps1_1, jyseps2_1, jyseps1_2, jyseps2_2,
   y_boundary_guards, exclude_y_guard_cells=False):
    ivertex_core = []
    Rxy_core = []
    Zxy_core = []
    if exclude_y_guard_cells:
        jyseps1_1g = jyseps1_1
        jyseps2_1g = jyseps2_1
    else:
        jyseps1_1g = jyseps1_1 + y_boundary_guards
        jyseps2_1g = jyseps2_1 + y_boundary_guards
    for j in range(jyseps1_1g+1,jyseps2_1g+1):
        ivertex_core.append(ivertex_corners[0,j])
        Rxy_core.append(data_Rxy[0,j])
        Zxy_core.append(data_Zxy[0,j])
    for j in range(jyseps1_2+1,jyseps2_2+1):
        ivertex_core.append(ivertex_corners[0,j])
        Rxy_core.append(data_Rxy[0,j])
        Zxy_core.append(data_Zxy[0,j])
    return ivertex_core, Rxy_core, Zxy_core

def get_target_ll_vertices(ivertex_corners, data_Rxy, data_Zxy,
                        ivertex_corners_lr, data_Rxy_lr, data_Zxy_lr):
    ivertex_target_ll = []
    Rxy_target_ll = []
    Zxy_target_ll = []
    Nx, Ny = ivertex_corners.shape
    j = 0
    for i in range(0,Nx):
        ivertex_target_ll.append(ivertex_corners[i,j])
        Rxy_target_ll.append(data_Rxy[i,j])
        Zxy_target_ll.append(data_Zxy[i,j])
    i = Nx - 1
    ivertex_target_ll.append(ivertex_corners_lr[i,j])
    Rxy_target_ll.append(data_Rxy_lr[i,j])
    Zxy_target_ll.append(data_Zxy_lr[i,j])
    return ivertex_target_ll, Rxy_target_ll, Zxy_target_ll

def get_target_ul_vertices(ivertex_corners_ul, data_Rxy_ul, data_Zxy_ul,
                        ivertex_corners_ur, data_Rxy_ur, data_Zxy_ur,
                        ny_inner):
    ivertex_target_ul = []
    Rxy_target_ul = []
    Zxy_target_ul = []
    Nx, Ny = ivertex_corners_ul.shape
    j = ny_inner - 1
    for i in range(0,Nx):
        ivertex_target_ul.append(ivertex_corners_ul[i,j])
        Rxy_target_ul.append(data_Rxy_ul[i,j])
        Zxy_target_ul.append(data_Zxy_ul[i,j])
    i = Nx - 1
    ivertex_target_ul.append(ivertex_corners_ur[i,j])
    Rxy_target_ul.append(data_Rxy_ur[i,j])
    Zxy_target_ul.append(data_Zxy_ur[i,j])
    return ivertex_target_ul, Rxy_target_ul, Zxy_target_ul

def get_target_ur_vertices(ivertex_corners, data_Rxy, data_Zxy,
                    ivertex_corners_lr, data_Rxy_lr, data_Zxy_lr,
                    ny_inner):
    ivertex_target_ur = []
    Rxy_target_ur = []
    Zxy_target_ur = []
    Nx, Ny = ivertex_corners.shape
    j = ny_inner
    for i in range(0,Nx):
        ivertex_target_ur.append(ivertex_corners[i,j])
        Rxy_target_ur.append(data_Rxy[i,j])
        Zxy_target_ur.append(data_Zxy[i,j])
    i = Nx - 1
    ivertex_target_ur.append(ivertex_corners_lr[i,j])
    Rxy_target_ur.append(data_Rxy_lr[i,j])
    Zxy_target_ur.append(data_Zxy_lr[i,j])
    return ivertex_target_ur, Rxy_target_ur, Zxy_target_ur

def get_target_lr_vertices(ivertex_corners_ur, data_Rxy_ur, data_Zxy_ur,
                        ivertex_corners_ul, data_Rxy_ul, data_Zxy_ul):
    ivertex_target_lr = []
    Rxy_target_lr = []
    Zxy_target_lr = []
    Nx, Ny = ivertex_corners_ur.shape
    j = Ny - 1
    for i in range(0,Nx):
        ivertex_target_lr.append(ivertex_corners_ur[i,j])
        Rxy_target_lr.append(data_Rxy_ur[i,j])
        Zxy_target_lr.append(data_Zxy_ur[i,j])
    i = 0
    ivertex_target_lr.append(ivertex_corners_ul[i,j])
    Rxy_target_lr.append(data_Rxy_ul[i,j])
    Zxy_target_lr.append(data_Zxy_ul[i,j])
    return ivertex_target_lr, Rxy_target_lr, Zxy_target_lr

def plot_corners_get_dmplex_data(file_path,interactive_plot=False,print_cells_to_screen_output=False):
    dataset = nc.Dataset(file_path)

    y_boundary_guards = np.copy(dataset.variables['y_boundary_guards'][...])
    ixseps1 = np.copy(dataset.variables['ixseps1'][...])
    ixseps2 = np.copy(dataset.variables['ixseps2'][...])
    jyseps1_1 = np.copy(dataset.variables['jyseps1_1'][...])
    jyseps2_1 = np.copy(dataset.variables['jyseps2_1'][...])
    jyseps1_2 = np.copy(dataset.variables['jyseps1_2'][...])
    jyseps2_2 = np.copy(dataset.variables['jyseps2_2'][...])
    ny_inner = np.copy(dataset.variables['ny_inner'][...])
    print("y_boundary_guards",y_boundary_guards)
    print("ixseps1",ixseps1)
    print("ixseps2",ixseps2)
    print("jyseps1_1",jyseps1_1)
    print("jyseps2_1",jyseps2_1)
    print("jyseps1_2",jyseps1_2)
    print("jyseps2_2",jyseps2_2)
    print("ny_inner",ny_inner)
    exclude_y_guard_cells = True
    if exclude_y_guard_cells:
        s = slice(0,-1,1), list(range(y_boundary_guards,ny_inner+y_boundary_guards)) + list(range(ny_inner+3*y_boundary_guards,2*ny_inner+3*y_boundary_guards))
    else:
        s = slice(0,-1,1), slice(0,-1,1)
    data_Rxy = np.copy(dataset.variables['Rxy_corners'][:])[s]
    data_Zxy = np.copy(dataset.variables['Zxy_corners'][:])[s]
    data_Rxy_lr = np.copy(dataset.variables['Rxy_lower_right_corners'][:])[s]
    data_Zxy_lr = np.copy(dataset.variables['Zxy_lower_right_corners'][:])[s]
    data_Rxy_ur = np.copy(dataset.variables['Rxy_upper_right_corners'][:])[s]
    data_Zxy_ur = np.copy(dataset.variables['Zxy_upper_right_corners'][:])[s]
    data_Rxy_ul = np.copy(dataset.variables['Rxy_upper_left_corners'][:])[s]
    data_Zxy_ul = np.copy(dataset.variables['Zxy_upper_left_corners'][:])[s]
    dataset.close()
    # check that the points are indeed unique
    unique_points_2D(data_Rxy,data_Zxy)
    unique_points_2D(data_Rxy_lr,data_Zxy_lr)
    unique_points_2D(data_Rxy_ur,data_Zxy_ur)
    unique_points_2D(data_Rxy_ul,data_Zxy_ul)
    # construct array with unique set of vertices covering the
    # Hypnotoad grid.
    Nx, Ny = data_Rxy.shape
    Npoint = Nx*Ny

    corners_Rxy = np.zeros((Nx+1,Ny+2))
    corners_Rxy[0:Nx,0:Ny//2] = data_Rxy[:,0:Ny//2]
    corners_Rxy[Nx,0:Ny//2] = data_Rxy_lr[-1,0:Ny//2]
    corners_Rxy[1:Nx+1,Ny//2] = data_Rxy_ur[:,Ny//2 - 1]
    corners_Rxy[0,Ny//2] = data_Rxy_ul[0,Ny//2 -1]

    corners_Rxy[0:Nx,Ny//2+1:Ny+1] = data_Rxy[:,Ny//2:Ny]
    corners_Rxy[Nx,Ny//2+1:Ny+1] = data_Rxy_lr[Nx-1,Ny//2:Ny]
    corners_Rxy[0:Nx,Ny+1] = data_Rxy_ul[0:Nx,Ny-1]
    corners_Rxy[Nx,Ny+1] = data_Rxy_ur[Nx-1,Ny-1]

    corners_Zxy = np.zeros((Nx+1,Ny+2))
    corners_Zxy[0:Nx,0:Ny//2] = data_Zxy[:,0:Ny//2]
    corners_Zxy[Nx,0:Ny//2] = data_Zxy_lr[-1,0:Ny//2]
    corners_Zxy[1:Nx+1,Ny//2] = data_Zxy_ur[:,Ny//2 - 1]
    corners_Zxy[0,Ny//2] = data_Zxy_ul[0,Ny//2 -1]

    corners_Zxy[0:Nx,Ny//2+1:Ny+1] = data_Zxy[:,Ny//2:Ny]
    corners_Zxy[Nx,Ny//2+1:Ny+1] = data_Zxy_lr[Nx-1,Ny//2:Ny]
    corners_Zxy[0:Nx,Ny+1] = data_Zxy_ul[0:Nx,Ny-1]
    corners_Zxy[Nx,Ny+1] = data_Zxy_ur[Nx-1,Ny-1]

    # check that the points are indeed unique
    unique_points_2D(corners_Rxy,corners_Zxy)

    # Visualise points with scatter plots
    # Get 1D lists of points for the scatter plot
    Rpoints = np.reshape(data_Rxy, (Npoint,))
    Zpoints = np.reshape(data_Zxy, (Npoint,))
    Rpoints_lr = np.reshape(data_Rxy_lr, (Npoint,))
    Zpoints_lr = np.reshape(data_Zxy_lr, (Npoint,))
    Rpoints_ur = np.reshape(data_Rxy_ur, (Npoint,))
    Zpoints_ur = np.reshape(data_Zxy_ur, (Npoint,))
    Rpoints_ul = np.reshape(data_Rxy_ul, (Npoint,))
    Zpoints_ul = np.reshape(data_Zxy_ul, (Npoint,))

    Npoint_full = (Nx+1)*(Ny+2)
    Rpoints_full = np.reshape(corners_Rxy, (Npoint_full,))
    Zpoints_full = np.reshape(corners_Zxy, (Npoint_full,))
    # remove any duplicate points
    Rpoints_full, Zpoints_full = remove_nonunique_points(Rpoints_full,Zpoints_full)
    unique_points_1D(Rpoints_full,Zpoints_full)

    save_ivertex_indices_to_netcdf(file_path,"test.nc",Rpoints_full,Zpoints_full)

    # get an array containing, for each cell, the list of vertices,
    # labelled by the global index defined implicitly by Rpoints_full, Zpoints_full
    # store this in `cell_vertices`, and store the RZ values of each vertex
    # cell_vertices[icell,ivertex] in R, Z = cell_vertices_RZ[icell,ivertex,:]
    cell_vertices, cell_vertices_RZ, vertex_list = get_cell_vertex_list(data_Rxy, data_Zxy,
            data_Rxy_lr, data_Zxy_lr, data_Rxy_ur, data_Zxy_ur,
            data_Rxy_ul, data_Zxy_ul, Rpoints_full, Zpoints_full)
    ncells, nvertices = np.shape(cell_vertices)
    if print_cells_to_screen_output:
        for icell in range(0,ncells):
            print(cell_vertices[icell,:])
            print(cell_vertices_RZ[icell,:,:])

    # test getting data for resaving to hypnotoad netcdf file
    ivertex_corners = get_boutxx_corner_index(data_Rxy,data_Zxy,Rpoints_full,Zpoints_full)
    ivertex_corners_lr = get_boutxx_corner_index(data_Rxy_lr,data_Zxy_lr,Rpoints_full,Zpoints_full)
    ivertex_corners_ul = get_boutxx_corner_index(data_Rxy_ul,data_Zxy_ul,Rpoints_full,Zpoints_full)
    ivertex_corners_ur = get_boutxx_corner_index(data_Rxy_ur,data_Zxy_ur,Rpoints_full,Zpoints_full)
    #print("ivertex_corners",ivertex_corners)
    #print("ivertex_corners_lr",ivertex_corners_lr)
    #print("ivertex_corners_ul",ivertex_corners_ul)
    #print("ivertex_corners_ur",ivertex_corners_ur)
    ivertex_pfr_lower, Rxy_pfr_lower, Zxy_pfr_lower = get_pfr_lower_boundary_vertices(ivertex_corners,data_Rxy,data_Zxy,
                                                        ivertex_corners_ul,data_Rxy_ul,data_Zxy_ul,
                                                        y_boundary_guards,jyseps1_1,exclude_y_guard_cells=exclude_y_guard_cells)

    ivertex_pfr_upper, Rxy_pfr_upper, Zxy_pfr_upper = get_pfr_upper_boundary_vertices(ivertex_corners,data_Rxy,data_Zxy,
                                                        ivertex_corners_ul,data_Rxy_ul,data_Zxy_ul,
                                                        y_boundary_guards,jyseps1_2,jyseps2_1,jyseps2_2,
                                                        ny_inner, exclude_y_guard_cells=exclude_y_guard_cells)

    ivertex_sol_vac_left, Rxy_sol_vac_left, Zxy_sol_vac_left = get_vac_left_boundary_vertices(ivertex_corners_lr,data_Rxy_lr,data_Zxy_lr,
                                                                    ivertex_corners_ur,data_Rxy_ur,data_Zxy_ur,
                                                                    y_boundary_guards,ny_inner)

    ivertex_sol_vac_right, Rxy_sol_vac_right, Zxy_sol_vac_right = get_vac_right_boundary_vertices(ivertex_corners_lr,data_Rxy_lr,data_Zxy_lr,
                                                                    ivertex_corners_ur,data_Rxy_ur,data_Zxy_ur,
                                                                    jyseps1_2,ny_inner)

    ivertex_core, Rxy_core, Zxy_core = get_core_boundary_vertices(ivertex_corners, data_Rxy, data_Zxy,
                                                                jyseps1_1, jyseps2_1, jyseps1_2, jyseps2_2,
                                                                y_boundary_guards, exclude_y_guard_cells=exclude_y_guard_cells)

    ivertex_target_ll, Rxy_target_ll, Zxy_target_ll = get_target_ll_vertices(ivertex_corners, data_Rxy, data_Zxy,
                                                                ivertex_corners_lr, data_Rxy_lr, data_Zxy_lr)

    ivertex_target_ul, Rxy_target_ul, Zxy_target_ul = get_target_ul_vertices(ivertex_corners_ul, data_Rxy_ul, data_Zxy_ul,
                                                                ivertex_corners_ur, data_Rxy_ur, data_Zxy_ur,
                                                                ny_inner)

    ivertex_target_ur, Rxy_target_ur, Zxy_target_ur = get_target_ur_vertices(ivertex_corners, data_Rxy, data_Zxy,
                                                                ivertex_corners_lr, data_Rxy_lr, data_Zxy_lr,
                                                                ny_inner)

    ivertex_target_lr, Rxy_target_lr, Zxy_target_lr = get_target_lr_vertices(ivertex_corners_ur, data_Rxy_ur, data_Zxy_ur,
                                                                ivertex_corners_ul, data_Rxy_ul, data_Zxy_ul)
    boundary_vertex_info = {
        "pfr_lower" : {"ivertex": ivertex_pfr_lower,
                       "Rxy" : Rxy_pfr_lower,
                       "Zxy" : Zxy_pfr_lower,
                       "color" : "g",
                       "marker" : "1",
                       "DMFaceSetsLabel" : 100,
                       },
        "pfr_upper" : {"ivertex": ivertex_pfr_upper,
                       "Rxy" : Rxy_pfr_upper,
                       "Zxy" : Zxy_pfr_upper,
                       "color" : "b",
                       "marker" : "2",
                       "DMFaceSetsLabel" : 100,
                       },
        "sol_vac_left" : {"ivertex": ivertex_sol_vac_left,
                       "Rxy" : Rxy_sol_vac_left,
                       "Zxy" : Zxy_sol_vac_left,
                       "color" : "r",
                       "marker" : "3",
                       "DMFaceSetsLabel" : 200,
                       },
        "sol_vac_right" : {"ivertex": ivertex_sol_vac_right,
                       "Rxy" : Rxy_sol_vac_right,
                       "Zxy" : Zxy_sol_vac_right,
                       "color" : "r",
                       "marker" : "3",
                       "DMFaceSetsLabel" : 200,
                       },
        "core" : {"ivertex": ivertex_core,
                       "Rxy" : Rxy_core,
                       "Zxy" : Zxy_core,
                       "color" : "k",
                       "marker" : "4",
                       "DMFaceSetsLabel" : 300,
                       },
        "target_ll" : {"ivertex": ivertex_target_ll,
                       "Rxy" : Rxy_target_ll,
                       "Zxy" : Zxy_target_ll,
                       "color" : "y",
                       "marker" : "3",
                       "DMFaceSetsLabel" : 400,
                       },
        "target_ul" : {"ivertex": ivertex_target_ul,
                       "Rxy" : Rxy_target_ul,
                       "Zxy" : Zxy_target_ul,
                       "color" : "y",
                       "marker" : "3",
                       "DMFaceSetsLabel" : 400,
                       },
        "target_ur" : {"ivertex": ivertex_target_ur,
                       "Rxy" : Rxy_target_ur,
                       "Zxy" : Zxy_target_ur,
                       "color" : "y",
                       "marker" : "3",
                       "DMFaceSetsLabel" : 400,
                       },
        "target_lr" : {"ivertex": ivertex_target_lr,
                       "Rxy" : Rxy_target_lr,
                       "Zxy" : Zxy_target_lr,
                       "color" : "y",
                       "marker" : "3",
                       "DMFaceSetsLabel" : 400,
                       },
    }

    # Make a scatter plot to show the mesh corners
    plt.figure(figsize=(10, 6))
    x = Rpoints
    y = Zpoints
    #scatter = plt.scatter(x, y, c='b',marker='1')
    #scatter = plt.scatter(Rpoints_lr, Zpoints_lr, c='r',marker='2')
    #scatter = plt.scatter(Rpoints_ur, Zpoints_ur, c='g',marker='3')
    #scatter = plt.scatter(Rpoints_ul, Zpoints_ul, c='k',marker='4')
    scatter = plt.scatter(Rpoints_full, Zpoints_full, c='m',marker='x')
    for key in boundary_vertex_info.keys():
        Rxy = boundary_vertex_info[key]["Rxy"]
        Zxy = boundary_vertex_info[key]["Zxy"]
        color = boundary_vertex_info[key]["color"]
        marker = boundary_vertex_info[key]["marker"]
        scatter = plt.scatter(Rxy,Zxy, c=color,marker=marker)
    # uncomment for labels on original data points
    #for ic in range(0,Npoint):
    #    plt.text(x[ic],y[ic],str(ic))
    #    plt.text(Rpoints_ul[ic],Zpoints_ul[ic],str(ic))
    # uncomment for labels on aggregated array of points
    #Npoint_full = len(Rpoints_full)
    #for ic in range(0,Npoint_full):
    #    plt.text(Rpoints_full[ic],Zpoints_full[ic],str(ic))
    for key in boundary_vertex_info.keys():
        Rxy = boundary_vertex_info[key]["Rxy"]
        Zxy = boundary_vertex_info[key]["Zxy"]
        ivertex = boundary_vertex_info[key]["ivertex"]
        Nvertex = len(ivertex)
        for i in range(0,Nvertex):
            plt.text(Rxy[i],Zxy[i],str(ivertex[i]))

    plt.title('Meshpoints')
    plt.xlabel('R')
    plt.ylabel('Z')
    if interactive_plot:
        plt.show()
    plt.savefig(file_path+".mesh_plot.pdf")
    return Nx, Ny, cell_vertices, vertex_list, boundary_vertex_info