import h5py
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from petsc4py import PETSc
import argparse
parser = argparse.ArgumentParser(description="Process a string input.")
parser.add_argument("dmplex_h5_file_path", type=str, help="The path to the HDF5 file representing the DMPlex data")
parser.add_argument("particle_trajectory_h5_file_path", type=str, help="The path to the HDF5 file representing the particle data")

args = parser.parse_args()
print(f"Animating particle paths from {args.particle_trajectory_h5_file_path} on DMPlex edges from {args.dmplex_h5_file_path}")

def load_dmplex(file_path):
    dm = PETSc.DMPlex().create()
    viewer = PETSc.Viewer().createHDF5(file_path, 'r')
    dm.load(viewer)
    viewer.destroy()
    dm.setFromOptions()
    dm.setUp()
    return dm

def get_mesh_edges(dm):
    # Get coordinates
    section = dm.getCoordinateSection()
    coords = dm.getCoordinatesLocal().array
    dim = dm.getCoordinateDim()
    # Get all edges (cone of each edge cell)
    edges = []
    for p in range(dm.getChart()[0], dm.getChart()[1]):
        # cone of edge will only have two points
        cone = dm.getCone(p)
        if len(cone) == 2:
            v0, v1 = cone
            x0 = coords[section.getOffset(v0):section.getOffset(v0)+dim]
            x1 = coords[section.getOffset(v1):section.getOffset(v1)+dim]
            edges.append((x0, x1))
    return edges

dm = load_dmplex(args.dmplex_h5_file_path)
#dm = load_dmplex('dmplex/expected_nonorthogonal.grd.nc.mesh.h5')
edges = get_mesh_edges(dm)

# Plot using object-oriented matplotlib
fig, ax = plt.subplots()

def load_particle_data(file_path):
    particle_data = h5py.File(file_path,"r")
    data_all_timesteps = list(particle_data.keys())
    nstep = len(data_all_timesteps)

    particle_positions = []
    for it in range(0,nstep):
        group = particle_data[f"Step#{it}"]
        #print(list(group.keys()))
        P_0 = group["P_0"]
        P_1 = group["P_1"]
        nparticles = len(P_0)
        pdata = np.zeros((nparticles,2))
        pdata[:,0] = P_0
        pdata[:,1] = P_1
        particle_positions.append(pdata)
    return particle_positions

def update_plot(i, data, scat):
    scat.set_offsets(data[i])
    return scat,

particle_positions = load_particle_data(args.particle_trajectory_h5_file_path)
nstep = len(particle_positions)

# plot an animation of the particles on the mesh
fig, ax = plt.subplots()
# Plot DMPlex edges
for x0, x1 in edges:
    ax.plot([x0[0], x1[0]], [x0[1], x1[1]], color='k', linewidth=0.5)
# Animate particles
scat = ax.scatter(particle_positions[0][:,0], particle_positions[0][:,1], c='b',s=1.0, marker='.')
ax.set_title('Particle Positions')
ax.set_xlabel('R')
ax.set_ylabel('Z')

def update(frame):
    scat.set_offsets(np.c_[particle_positions[frame][:,0], particle_positions[frame][:,1]])
    return scat,

ani = FuncAnimation(fig, update, frames=nstep, interval=50, blit=True)
output_path = args.particle_trajectory_h5_file_path + ".animation.gif"
ani.save(output_path)
print(f"Saving animation of particle paths to {output_path}")
#plt.show()




