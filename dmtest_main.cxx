// #include "bout/mpi_wrapper.hxx"
#include "bout/petsclib.hxx"
#include <petscdmplex.h>
#include <petscsys.h>
#include <petscviewer.h>
#include <vector>
#include <cmath>
// N.B. comment the remaining includes to permit
// compilation without errors
#include "bout/field_factory.hxx"
#include "bout/bout.hxx"
#include "bout/version.hxx"
#include "div_ops.hxx"
#include "bout/fv_ops.hxx"
#include "bout/difops.hxx"
#include "../include/div_ops.hxx"

int main(int argc, char **argv) {
    // attempt to call BOUT to 
    // get information from a BOUT
    // mesh object
    // N.B. Comment the next three lines
    // to permit compilation as is
    BoutInitialise(argc, argv);
    Mesh* mesh = Mesh::create(&Options::root()["mesh"]);
    mesh->load();

    // Now a program to make a custom DMPlex object
    // and prepare to visualise it with a .vtk output file
    // Based on 
    // https://github.com/will-saunders-ukaea/petsc_neso_particles_example/tree/main#
    PetscInitialize(&argc, &argv, NULL, NULL);
    double pi = std::acos(-1.0);
    // First we setup the topology of the mesh.
    PetscInt num_cells_per_side = 50;
    PetscInt num_cells_owned = num_cells_per_side * num_cells_per_side;
    // PetscInt mpi_rank = 0;
    // PetscInt mpi_size = num_cells_owned;
    std::vector<PetscInt> cells;
    cells.reserve(num_cells_owned * 4);

    // We are careful to list the vertices in counter clock-wise order, this might
    // matter.
    for (PetscInt cy = 0; cy < num_cells_per_side; cy++) {
        for (PetscInt cx = 0; cx < num_cells_per_side; cx++) {
            // These are global indices not local indices.
            PetscInt vertex_sw = cx + cy * (num_cells_per_side + 1);
            PetscInt vertex_se = vertex_sw + 1;
            PetscInt vertex_ne = vertex_se + num_cells_per_side + 1;
            PetscInt vertex_nw = vertex_ne - 1;
            cells.push_back(vertex_sw);
            cells.push_back(vertex_se);
            cells.push_back(vertex_ne);
            cells.push_back(vertex_nw);
        }
    }
    /*
    * Each rank owns a contiguous block of global indices. We label our indices
    * lexicographically (row-wise). Sorting out the global vertex indexing is
    * probably one of the more tedious parts.
    */
    PetscInt num_vertices_owned = (num_cells_per_side + 1)*(num_cells_per_side + 1);
    
    /*
    * Create the coordinates for the block of vertices we pass to petsc. For an
    * existing mesh in memory this step will probably involve some MPI
    * communication to gather the blocks of coordinates on the ranks which pass
    * them to PETSc.
    */
    double twopi_N = 2.0*pi/num_cells_per_side;
    PetscInt pc;
    std::vector<PetscScalar> vertex_coords(num_vertices_owned * 2);
    for (PetscInt py = 0; py < num_cells_per_side + 1; py++) {
        for (PetscInt px = 0; px < num_cells_per_side + 1; px++) {
            // a compound index to index the vertices
            pc = px + py * (num_cells_per_side + 1);
            // our cell extent is 1.0
            vertex_coords.at(pc * 2 + 0) = (px+1)*std::cos(twopi_N*py);
            vertex_coords.at(pc * 2 + 1) = (px+1)*std::sin(twopi_N*py);
        }
    }
    // This DM will contain the DMPlex after we call the creation routine.
    DM dm;
    // Create the DMPlex from the cells and coordinates.
    PetscErrorCode ierr = DMPlexCreateFromCellListParallelPetsc(
        PETSC_COMM_WORLD, 2, num_cells_owned, num_vertices_owned, PETSC_DECIDE, 4,
        PETSC_TRUE, cells.data(), 2, vertex_coords.data(), NULL, NULL, &dm);
    CHKERRQ(ierr);
  
    // Distribute the mesh if running in parallel
    DM distributedDM = NULL;
    ierr = DMPlexDistribute(dm, 0, NULL, &distributedDM); CHKERRQ(ierr);
    if (distributedDM) {
        DMDestroy(&dm);
        dm = distributedDM;
    }

    // View the mesh in VTK format
    PetscViewer viewer;
    ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD, "mesh.vtk", FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
    ierr = DMView(dm, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    // Clean up
    DMDestroy(&dm);
    PetscFinalize();
    return 0;
}
