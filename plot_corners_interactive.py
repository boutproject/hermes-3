from plot_corners_functions import plot_corners_get_dmplex_data

file_path = 'dmtest_data/expected_nonorthogonal.grd.nc'
Nx,Ny,cell_vertices,vertex_list = plot_corners_get_dmplex_data(file_path,interactive_plot=True,print_cells_to_screen_output=True)