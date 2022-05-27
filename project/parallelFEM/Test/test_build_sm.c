#include "hpc.h"
#include <math.h>

#ifndef MPI_GRID_X 
#define MPI_GRID_X  2
#endif

#ifndef MPI_GRID_Y
#define MPI_GRID_Y  2
#endif

#ifndef NUM_MESH_REFINES
#define NUM_MESH_REFINES 2
#endif

double F_vol(double x[2], index typ) { return (0.0); }

double g_Neu(double x[2], index typ) { return (x[0] * x[1]); }

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
	printf("\n=== Start test_mesh_transfer ==\n");	
	
	// Load basic mesh
	mesh *m1 = mesh_load(fname);
	mesh_getEdge2no(m1->nelem, m1->elem, &m1->nedges, &m1->edge2no);
	m1->fixed = mesh_getFixed(m1->ncoord, m1->bdry, m1->nbdry, &m1->nfixed);
	
	// initial refinement
	printf("Refining mesh\n");
	meshRefined = mesh_initRefinement(m1, NUM_MESH_REFINES); 
	printf("New mesh dimension (%td,%td)\n",(index)(sqrt((double)meshRefined->ncoord)),
		(index)(sqrt((double)meshRefined->ncoord)));
	
	// create mapping
	mapping = mesh_split(meshRefined, dims);

	// transfer local meshes
	localMapping = mesh_transfer(mapping, grid);

    // Nonrot processes: Call method mesh_transfer to receive mesh and mapping
    // data
    } else{
	localMapping = mesh_transfer(mapping, grid);
    }

    
    // Build the local stiffness matrices
    printf("Rank %d building local stiffness matrix\n",rank);
    mesh *localMesh = localMapping->localMesh;
    mesh_getEdge2no(localMesh->nelem, localMesh->elem, &localMesh->nedges, &localMesh->edge2no);
    localMesh->fixed = mesh_getFixed(localMesh->ncoord, localMesh->bdry, localMesh->nbdry, &localMesh->nfixed);
    sed *sm_local = sed_sm_build(localMesh);

    // Display local stiffness matrix if rank = 1
    if (rank==1){
	sed_print(sm_local,0);
    }

    // build local rhs
    double *b = malloc(localMesh->ncoord*sizeof(double));
    mesh_build_rhs(localMesh, b, F_vol, g_Neu); 






    // If rank=0 free the mapping, the refined mesh and the base mesh
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

    if (rank==0){
	printf("\n=== End test_mesh_transfer ===\n");
    }

    MPI_Finalize();
}
