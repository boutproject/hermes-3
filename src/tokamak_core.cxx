#include "../include/tokamak_core.hxx"

#include <bout/mesh.hxx>
using bout::globals::mesh;

void TokamakCore::transform(Options& state) {
  Options& species = state["species"][name];

  Coordinates* coord = mesh->getCoordinates();
  const Field2D& J = coord->J;
  const Field2D& dy = coord->dy;
  const Field2D& dx = coord->dx;
  const Field2D& dz = coord->dz;

  energy_source = 0;
  density_source = 0;

  /// Assign core sources
  if ((power > 0) or (particle_flow > 0)) {

    core_volume = 0;
    core_volume_local = 0;
    
    /// Loop through all core cells (excl.guards) and tally up their volume for this processor
    if(mesh->firstX()){   // Only do this for the processor which has the core region
      if (mesh->periodicY(mesh->xstart)) {   // Only do this for the processor with a periodic Y (core)
        for(int iy=mesh->ystart; iy<=mesh->yend; iy++){
          for(int iz=mesh->zstart; iz<=mesh->zend; iz++){

            core_volume_local += J(mesh->xstart, iy) * dx(mesh->xstart, iy) * dy(mesh->xstart, iy) * dz(mesh->xstart, iy);

          }
        }
      }
    }

    // Sum over all processors 
    // https://mpitutorial.com/tutorials/mpi-reduce-and-allreduce/
    MPI_Allreduce(&core_volume_local, &core_volume, 1, MPI_DOUBLE, MPI_SUM, BoutComm::get());

    /// Calculate the source and assign it to all cells in the core ring
    if(mesh->firstX()){ 
      if (mesh->periodicY(mesh->xstart)) {  
        for(int iy=0; iy < mesh->LocalNy ; iy++){
          for(int iz=0; iz < mesh->LocalNz; iz++){

              energy_source(mesh->xstart, iy, iz) = power / core_volume;
              density_source(mesh->xstart, iy, iz) = particle_flow / core_volume;
            
          }
        }
      }
    }

    // Put the sources in the state so that evolve_pressure picks them up
    set<Field3D>(species["energy_source"], energy_source);
    set<Field3D>(species["density_source"], density_source);
  }
}
