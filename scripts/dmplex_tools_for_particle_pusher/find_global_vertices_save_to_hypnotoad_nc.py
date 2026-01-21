from plot_corners_functions import plot_corners_get_dmplex_data
import argparse
parser = argparse.ArgumentParser(description="Process a Hypnotoad mesh file to produce a modified mesh file and a file containing the global list of vertices of the mesh. Produce plots of cell vertices.")
parser.add_argument("hypnotoad_nc_file_path", type=str, help="The path to the Hypnotoad netCDF file representing the mesh.")
args = parser.parse_args()
file_path = args.hypnotoad_nc_file_path
print(f"Processing Hypnotoad mesh from {args.hypnotoad_nc_file_path}")

Nx,Ny,cell_vertices,vertex_list,boundary_vertex_info = plot_corners_get_dmplex_data(file_path,interactive_plot=True,print_cells_to_screen_output=False)
