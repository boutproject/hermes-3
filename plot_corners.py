import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

def isapprox(a,b,tol=1.0e-12):
        if abs(a - b) < tol:
            return True
        else:
            return False

def unique_points(corners_Rxy,corners_Zxy):
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
unique_points(data_Rxy,data_Zxy)
unique_points(data_Rxy_lr,data_Zxy_lr)
unique_points(data_Rxy_ur,data_Zxy_ur)
unique_points(data_Rxy_ul,data_Zxy_ul)
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
unique_points(corners_Rxy,corners_Zxy)

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

x = Rpoints
y = Zpoints
# Make a scatter plot to show the mesh corners
plt.figure(figsize=(10, 6))
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
for ic in range(0,Npoint_full):
    plt.text(Rpoints_full[ic],Zpoints_full[ic],str(ic))
plt.title('Meshpoints')
plt.xlabel('R')
plt.ylabel('Z')
plt.show()
