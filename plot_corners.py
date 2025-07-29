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
# Get 1D lists of points for the scatter plot
Rpoints = np.reshape(data_Rxy, (Npoint,))
Zpoints = np.reshape(data_Zxy, (Npoint,))
Rpoints_lr = np.reshape(data_Rxy_lr, (Npoint,))
Zpoints_lr = np.reshape(data_Zxy_lr, (Npoint,))
Rpoints_ur = np.reshape(data_Rxy_ur, (Npoint,))
Zpoints_ur = np.reshape(data_Zxy_ur, (Npoint,))
Rpoints_ul = np.reshape(data_Rxy_ul, (Npoint,))
Zpoints_ul = np.reshape(data_Zxy_ul, (Npoint,))
x = Rpoints
y = Zpoints
# Make a scatter plot to show the mesh corners
plt.figure(figsize=(10, 6))
scatter = plt.scatter(x, y, c='b',marker='1')
scatter = plt.scatter(Rpoints_lr, Zpoints_lr, c='r',marker='2')
scatter = plt.scatter(Rpoints_ur, Zpoints_ur, c='g',marker='3')
scatter = plt.scatter(Rpoints_ul, Zpoints_ul, c='k',marker='4')
plt.title('Meshpoints')
plt.xlabel('R')
plt.ylabel('Z')
plt.show()
