from plot_corners_functions import plot_corners_get_dmplex_data

file_path = 'dmtest_data/expected_nonorthogonal.grd.nc'
#file_path = 'dmtest_data/nonorthog.bout.grd.nc'
#file_path = 'dmtest_data/lsn.bout.grd.nc'
Nx,Ny,cell_vertices,vertex_list,boundary_vertex_info = plot_corners_get_dmplex_data(file_path,interactive_plot=True,print_cells_to_screen_output=False)
