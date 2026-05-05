import zoidberg 
from zoidberg.field import Slab
import numpy as np

import os
import sys

class Axial_circular(Slab):
    def __init__(self,rho_1,rho_2,Lz,By,R0, q0 = 3.0, q1 = 0.2):
        self.rho_1 = float(rho_1)
        self.rho_2 = float(rho_2)
        self.Lz = float(Lz)
        self.By = float(By)
        self.R0 = R0
        
    def qprofile(self,x,z,phi):
        theta=np.arctan2(z,x-self.R0)
        r = np.sqrt(np.power(x-self.R0,2)+np.power(z,2))
        return q0 + q1 * r
        
    def Bxfunc(self,x,z,phi):
        theta=np.arctan2(z,x-self.R0)
        r = np.sqrt(np.power(x-self.R0,2)+np.power(z,2))
        q = self.qprofile(x,z,phi)
        return -r * np.sin(theta) / q
    
    def Bzfunc(self,x,z,phi):
        theta=np.arctan2(z,x-self.R0)
        r = np.sqrt(np.power(x-self.R0,2)+np.power(z,2))
        q = self.qprofile(x,z,phi)
        return r * np.cos(theta) / q   

    def Byfunc(self,x,z,phi):
        return np.full(x.shape,self.By)




R0 = 1.5
rho_1 = 0.4
rho_2 = 0.6
B0 = 1.0
Ly = 2.0 * np.pi
q0 = 3.0
q1 = 0.0
field = Axial_circular(rho_1, rho_2, Ly, B0,R0, q0=q0, q1=q1)

folder = "mms_axial_circular_fci/"

if not os.path.exists(folder):
    os.mkdir(folder)

force = "-f" in sys.argv or "--force" in sys.argv


    
for scale in [1,2]:
    nx = scale * 16
    nz = scale * 64
    ny = scale * 8

    filename = f"Axial_circular_q_{nx}_{nz}_{rho_1}_{rho_2}_{q0}_{q1}.fci.grid.nc"
    fn = folder + filename
    
    if os.path.exists(fn) and not force:
        print(fn, " exists")
        continue
    
    
    inner = zoidberg.rzline.shaped_line(R0=R0,a=rho_1,n=nz)
    outer = zoidberg.rzline.shaped_line(R0=R0,a=rho_2,n=nz)

    pol_grid = zoidberg.poloidal_grid.grid_annulus(inner,outer,nx+4,nz)
    pol_grids = []
    for i in range(ny):
        pol_grids.append(pol_grid)
    
    grid = zoidberg.grid.Grid(pol_grids, np.linspace(0.0, 2.0 * np.pi, ny, endpoint=False),Ly,yperiodic=True)
    maps = zoidberg.zoidberg.make_maps(grid, field)

    maps["forward_xt_prime"][2,:,:] = 2.0
    maps["backward_xt_prime"][2,:,:] = 2.0

    maps["forward_xt_prime"][-3,:,:] = float(nx+4-3)
    maps["backward_xt_prime"][-3,:,:] = float(nx+4-3)

    with zoidberg.zoidberg.MapWriter(fn) as mw:
        mw.add_grid_field(grid, field)
        mw.add_maps(maps)
        mw.add_dagp()

    print(fn, " done")
