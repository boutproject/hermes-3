import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

def isapprox(a,b,tol=1.0e-12):
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
def convert_R_Z_to_vertex_index(R,Z,R_vertices,Z_vertices):
    Nc = len(R_vertices)
    for ic in range(0,Nc):
        if isapprox(R,R_vertices[ic]) and isapprox(Z,Z_vertices[ic]):
            ic_vertex = ic
            return ic_vertex
    print("Failed to find (R,Z) index")
    return None

def get_cell_vertex_list_inner(R_ll,Z_ll, R_lr, Z_lr,
                        R_ur, Z_ur, R_ul, Z_ul,
                        R_vertices, Z_vertices):
    vertices = []
    # go in anti-clockwise order
    vertices.append(convert_R_Z_to_vertex_index(R_ll,Z_ll,R_vertices,Z_vertices))
    vertices.append(convert_R_Z_to_vertex_index(R_lr,Z_lr,R_vertices,Z_vertices))
    vertices.append(convert_R_Z_to_vertex_index(R_ur,Z_ur,R_vertices,Z_vertices))
    vertices.append(convert_R_Z_to_vertex_index(R_ul,Z_ul,R_vertices,Z_vertices))
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
    return cell_vertices, cell_vertices_RZ


file_path = 'dmtest_data/expected_nonorthogonal.grd.nc'
#file_path = 'dmtest_data/expected_orthogonal.grd.nc'
dataset = nc.Dataset(file_path)

data_Rxy = dataset.variables['Rxy_corners'][:]
data_Zxy = dataset.variables['Zxy_corners'][:]
data_Rxy_lr = dataset.variables['Rxy_lower_right_corners'][:]
data_Zxy_lr = dataset.variables['Zxy_lower_right_corners'][:]
data_Rxy_ur = dataset.variables['Rxy_upper_right_corners'][:]
data_Zxy_ur = dataset.variables['Zxy_upper_right_corners'][:]
data_Rxy_ul = dataset.variables['Rxy_upper_left_corners'][:]
data_Zxy_ul = dataset.variables['Zxy_upper_left_corners'][:]

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

# get an array containing, for each cell, the list of vertices,
# labelled by the global index defined implicitly by Rpoints_full, Zpoints_full
# store this in `cell_vertices`, and store the RZ values of each vertex
# cell_vertices[icell,ivertex] in R, Z = cell_vertices_RZ[icell,ivertex,:]
cell_vertices, cell_vertices_RZ = get_cell_vertex_list(data_Rxy, data_Zxy,
        data_Rxy_lr, data_Zxy_lr, data_Rxy_ur, data_Zxy_ur,
        data_Rxy_ul, data_Zxy_ul, Rpoints_full, Zpoints_full)
ncells, nvertices = np.shape(cell_vertices)
for icell in range(0,ncells):
    print(cell_vertices[icell,:])
    print(cell_vertices_RZ[icell,:,:])


# Make a scatter plot to show the mesh corners
plt.figure(figsize=(10, 6))
x = Rpoints
y = Zpoints
scatter = plt.scatter(x, y, c='b',marker='1')
scatter = plt.scatter(Rpoints_lr, Zpoints_lr, c='r',marker='2')
scatter = plt.scatter(Rpoints_ur, Zpoints_ur, c='g',marker='3')
scatter = plt.scatter(Rpoints_ul, Zpoints_ul, c='k',marker='4')
scatter = plt.scatter(Rpoints_full, Zpoints_full, c='m',marker='x')
# uncomment for labels on original data points
#for ic in range(0,Npoint):
#    plt.text(x[ic],y[ic],str(ic))
#    plt.text(Rpoints_ul[ic],Zpoints_ul[ic],str(ic))
# uncomment for labels on aggregated array of points
Npoint_full = len(Rpoints_full)
for ic in range(0,Npoint_full):
    plt.text(Rpoints_full[ic],Zpoints_full[ic],str(ic))
plt.title('Meshpoints')
plt.xlabel('R')
plt.ylabel('Z')
plt.show()
