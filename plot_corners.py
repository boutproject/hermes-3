import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

file_path = 'dmtest_data/expected_orthogonal.grd.nc'
dataset = nc.Dataset(file_path)

data_Rxy = dataset.variables['Rxy_corners'][:]
data_Zxy = dataset.variables['Zxy_corners'][:]
data_Rxy_lr = dataset.variables['Rxy_lower_right_corners'][:]
data_Zxy_lr = dataset.variables['Zxy_lower_right_corners'][:]
data_Rxy_ur = dataset.variables['Rxy_upper_right_corners'][:]
data_Zxy_ur = dataset.variables['Zxy_upper_right_corners'][:]
data_Rxy_ul = dataset.variables['Rxy_upper_left_corners'][:]
data_Zxy_ul = dataset.variables['Zxy_upper_left_corners'][:]

Nx, Ny = data_Rxy.shape
Npoint = Nx*Ny

corners_Rxy = np.zeros((Nx+1,Ny+1))
corners_Rxy[0:Nx,0:Ny] = data_Rxy[:,:]
corners_Rxy[Nx,0:Ny] = data_Rxy_lr[-1,:]
corners_Rxy[0:Nx,Ny] = data_Rxy_ur[:,-1]
corners_Rxy[Nx,Ny] = data_Rxy_ul[-1,-1]
corners_Zxy = np.zeros((Nx+1,Ny+1))
corners_Zxy[0:Nx,0:Ny] = data_Zxy[:,:]
corners_Zxy[Nx,0:Ny] = data_Zxy_lr[-1,:]
corners_Zxy[0:Nx,Ny] = data_Zxy_ur[:,-1]
corners_Zxy[Nx,Ny] = data_Zxy_ul[-1,-1]
#corners_Rxy[0:Nx,Ny] = data_Rxy_
# Get 1D lists of points for the scatter plot
Rpoints = np.reshape(data_Rxy, (Npoint,))
Zpoints = np.reshape(data_Zxy, (Npoint,))
Rpoints_lr = np.reshape(data_Rxy_lr, (Npoint,))
Zpoints_lr = np.reshape(data_Zxy_lr, (Npoint,))
Rpoints_ur = np.reshape(data_Rxy_ur, (Npoint,))
Zpoints_ur = np.reshape(data_Zxy_ur, (Npoint,))
Rpoints_ul = np.reshape(data_Rxy_ul, (Npoint,))
Zpoints_ul = np.reshape(data_Zxy_ul, (Npoint,))

Npoint_full = (Nx+1)*(Ny+1)
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
#for ic in range(0,Npoint):
#    plt.text(x[ic],y[ic],str(ic))
#    plt.text(Rpoints_lr[ic],Zpoints_lr[ic],str(ic))
for ic in range(0,Npoint_full):
    plt.text(Rpoints_full[ic],Zpoints_full[ic],str(ic))
plt.title('Meshpoints')
plt.xlabel('R')
plt.ylabel('Z')
plt.show()
