#include <iostream>
#include <vector>
#include <netcdf>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

int main() {
    try {
        // Open the NetCDF file
        NcFile dataFile("expected_orthogonal.grd.nc", NcFile::read);

        // Get the variables at the sw, se, ne, nw corners of cells
        NcVar Zxy_sw = dataFile.getVar("Zxy_corners");
        NcVar Zxy_se = dataFile.getVar("Zxy_lower_right_corners");
        NcVar Zxy_ne = dataFile.getVar("Zxy_upper_right_corners");
        NcVar Zxy_nw = dataFile.getVar("Zxy_upper_left_corners");
        NcVar Rxy_sw = dataFile.getVar("Rxy_corners");
        NcVar Rxy_se = dataFile.getVar("Rxy_lower_right_corners");
        NcVar Rxy_ne = dataFile.getVar("Rxy_upper_right_corners");
        NcVar Rxy_nw = dataFile.getVar("Rxy_upper_left_corners");

        // Get the dimensions
        vector<NcDim> dims = Rxy_sw.getDims();
        int Nx = dims[0].getSize();
        int Ny = dims[1].getSize();

        std::vector<float> Rxy_sw_data(Nx*Ny);
        std::vector<float> Rxy_se_data(Nx*Ny);
        std::vector<float> Rxy_ne_data(Nx*Ny);
        std::vector<float> Rxy_nw_data(Nx*Ny);
        std::vector<float> Zxy_sw_data(Nx*Ny);
        std::vector<float> Zxy_se_data(Nx*Ny);
        std::vector<float> Zxy_ne_data(Nx*Ny);
        std::vector<float> Zxy_nw_data(Nx*Ny);
        // load the data into buffers
        Rxy_sw.getVar(Rxy_sw_data.data());
        Rxy_se.getVar(Rxy_se_data.data());
        Rxy_ne.getVar(Rxy_ne_data.data());
        Rxy_nw.getVar(Rxy_nw_data.data());
        Zxy_sw.getVar(Zxy_sw_data.data());
        Zxy_se.getVar(Zxy_se_data.data());
        Zxy_ne.getVar(Zxy_ne_data.data());
        Zxy_nw.getVar(Zxy_nw_data.data());

        cout << "Nx = " << Nx << endl;
        cout << "Ny = " << Ny << endl;
        int ic = 0;
        for (int ix = 0; ix < Nx; ix++){
            for (int iy = 0; iy < Ny; iy++){
                ic = iy + Ny*ix;
                cout << " (iy,ix,ic) = " << iy << " " << ix << " " << ic << endl;
                cout << " Rxy(sw,se,ne,nw) = " << Rxy_sw_data[ic] << " " << Rxy_se_data[ic] << " " << Rxy_ne_data[ic] << " " << Rxy_nw_data[ic] << endl;
                cout << " Zxy(sw,se,ne,nw) = " << Zxy_sw_data[ic] << " " << Zxy_se_data[ic] << " " << Zxy_ne_data[ic] << " " << Zxy_nw_data[ic] << endl;
            }
        }
        
    } catch (NcException& e) {
        cerr << "NetCDF error: " << e.what() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
