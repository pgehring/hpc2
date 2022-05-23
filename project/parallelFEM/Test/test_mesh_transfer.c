#include <mpi.h>
#include "hpc.h"

#ifndef MPI_GRID_X 
#define MPI_GRID_X  2
#endif

#ifndef MPI_GRID_Y
#define MPI_GRID_Y  2
#endif

int main(int argc, char**argv){
    MPI_Init(&argc, &argv);

    // Initialize MPI gird Communication
    int nof_processes; MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int dims[2] = {MPI_GRID_X,MPI_GRID_Y};
    int periods[2] = {false,false};
    
    bool resultTotal = 1;

    char fname[32] = "../Problem/problem1";

    MPI_Dims_create(nof_processes, 2, dims);
    MPI_Comm grid;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &grid);
    MPI_Comm_rank(grid, &rank);


    MeshMapping ***mapping;
    mesh *m1;
    mesh *meshRefined;
    MeshMapping *localMapping;
    
    // root: load and refine mesh, send local meshes
    if (rank==0){
	DEBUG_PRINT("\n=== Start test_mesh_transfer ==\n");	
	
	// Load basic mesh
	mesh *m1 = mesh_load(fname);
	mesh_getEdge2no(m1->nelem, m1->elem, &m1->nedges, &m1->edge2no);
	m1->fixed = mesh_getFixed(m1->ncoord, m1->bdry, m1->nbdry, &m1->nfixed);
	
	// initial refinement
	meshRefined = mesh_initRefinement(m1, HPC_MAX(dims[0],dims[1]));

	// create mapping
	mapping = mesh_split(meshRefined, dims);

	// transfer local meshes
	localMapping = mesh_transfer(mapping, grid);
    } else{
	localMapping = mesh_transfer(mapping, grid);
    }

    // If rank=0 free the mapping the refined mesh and the base mesh
    if (rank ==0){
	delete2DMeshMapping(mapping, dims[0]);
	mesh_free(meshRefined);
	mesh_free(m1);

    // if rank!=0 free the local mesh
    } else{
	if (localMapping != NULL){
	    mesh_free(localMapping->localMesh);
	    deleteMeshMapping(localMapping);
	}
    }

    DEBUG_PRINT("\nProcess with rank %d ending Test!\n",rank);    

    MPI_Finalize();
}
