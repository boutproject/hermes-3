from petsc4py import PETSc
import argparse
parser = argparse.ArgumentParser(description="Process a string input.")
parser.add_argument("dmplex_h5_file_path", type=str, help="The path to the HDF5 file representing the DMPlex data")

args = parser.parse_args()
print(f"Plotting: {args.dmplex_h5_file_path}")
file_path = args.dmplex_h5_file_path

dm = PETSc.DMPlex().create()
viewer = PETSc.Viewer().createHDF5(file_path, 'r')
dm.load(viewer)
viewer.destroy()
dm.setFromOptions()
dm.setUp()

import matplotlib.pyplot as plt


# Get coordinates
section = dm.getCoordinateSection()
coords = dm.getCoordinatesLocal().array
dim = dm.getCoordinateDim()

# Get all edges (cone of each edge cell)
edges = []
for p in range(dm.getChart()[0], dm.getChart()[1]):
    #if dm.getLabelValue("celltype", p) == PETSc.DMPlex.CellType.EDGE:
    cone = dm.getCone(p)
    if len(cone) == 2:
        v0, v1 = cone
        x0 = coords[section.getOffset(v0):section.getOffset(v0)+dim]
        x1 = coords[section.getOffset(v1):section.getOffset(v1)+dim]
        edges.append((x0, x1))

# Plot using object-oriented matplotlib
fig, ax = plt.subplots()
for x0, x1 in edges:
    ax.plot([x0[0], x1[0]], [x0[1], x1[1]], color='k', linewidth=0.5)

ax.set_aspect('equal')
ax.set_title('DMPlex Mesh Edges')
ax.set_xlabel('x')
ax.set_ylabel('y')
output_path = file_path[:-3]+'.pdf'
plt.savefig(output_path)
print(f"Saving file to: {output_path}")
