#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hpc.h"


#define TIME_SAVE(j,tv) (gettimeofday(&tv[j], (struct timezone *)0))
#define TIME_ELAPSED(j, k, tv) \
    (1.E+6 * (tv[k].tv_sec - tv[j].tv_sec) + (tv[k].tv_usec - tv[j].tv_usec))

double kappa(double x[2], index typ)
{
    return (1.0);
}

double F_vol(double x[2], index typ) { return (0.0); }

double g_Neu(double x[2], index typ) { return (x[0] * x[1]); }

double u_D(double x[2])
{
    //  return ( 0.0 );
    return (x[0] * x[1]);
}

void initResultFile(char * fname, FILE *fileResult){
    fileResult = fopen(fname,"r");
   
    printf("Testing if file exists or is empty\n");
    // write header if result file is empty
    if (fileResult == NULL || fgetc(fileResult) == EOF){
	fclose(fileResult);
	fileResult = fopen(fname, "w+");
        fprintf(fileResult, "%10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
		 "nof_proc", "DOFs","Iterations", "T_Mesh","T_Mapping","T_Transfer",
		 "T_Setup","T_Solv","T_Res");
	fclose(fileResult);
    }else {
	fclose(fileResult);
    }
}


void saveBenchResults(struct timeval *tv, int n, index nofDoF, int nofIt, 
		      int nof_processes, FILE *benchFile){
    double durationVals[n-1];
    for (int i=1; i<n; ++i){
	durationVals[i-1] = 1.E+3*(tv[i].tv_sec - tv[0].tv_sec)+1.E-3*(tv[i].tv_usec -
		       tv[0].tv_usec);
    }
    
    fprintf(benchFile,"%10d %10td %10td ",nof_processes, nofDoF, nofIt); 

    for (int i=0; i<n-2; ++i){
	printf("T%d: %lf\n",i+1, durationVals[i]);
	fprintf(benchFile, "%10lf ",durationVals[i]);
    }
    
    printf("T%d: %lf\n", n-1, durationVals[n-2]);
    fprintf(benchFile, "%10lf\n",n-2, durationVals[n-2]);
}


/** Wrapper function for solving the demo problem with test implementation
    of CG solver */
int benchPoissonCG(char *fname, int numRefines, 
		       double (*fV)(double *, index), double (*fN)(double *, index),
		       struct timeval *tv, MPI_Comm grid){

    // Get MPI Rank and dimensions of grid
    int rank, nof_processes, gridDims[2], periods[2], coords[2];
    MPI_Comm_rank(grid, &rank);
    MPI_Comm_size(grid, &nof_processes);
    MPI_Cart_get(grid, 2, gridDims, periods, coords);

    // Declare pointers for mapping structs
    MeshMapping ***mapping;
    mesh *m1;
    mesh *meshRefined;
    MeshMapping *localMapping;

    // Declare pointer to solution vectors and right hand side
    double *localSolCG, *glblSolCG, *rhs;

    // number of iterations of solver
    int nofIt;

    // Save initial time value
    if (rank==0) TIME_SAVE(0,tv);


    // root: load and refine mesh, send local meshes
    if (rank==0){
	
	mesh *m1 = mesh_load(fname);
	mesh_getEdge2no(m1->nelem, m1->elem, &m1->nedges, &m1->edge2no);
	
	meshRefined = mesh_initRefinement(m1, numRefines); 
	// Save time after mesh loading and refinement
	TIME_SAVE(1,tv);
	
	mapping = mesh_split(meshRefined, gridDims);

	// Save time after mesh split
	TIME_SAVE(2,tv);
	localMapping = mesh_transfer(mapping, grid);

	/* Save time after all process got their loacl mesh (possible because
	    of blocking communication */
	TIME_SAVE(3, tv);

    } else{
	localMapping = mesh_transfer(mapping, grid);
    }
    
    // problem setup (stiffness matrix rhs and dirichlet data)
    mesh_getFixedNodes(localMapping);
    
    mesh *localMesh = localMapping->localMesh;
    sed *sm_local = sed_sm_build(localMesh);
   
    rhs = newVectorWithInit(localMesh->ncoord);
    mesh_build_rhs(localMesh, rhs, F_vol, g_Neu);
    
    localSolCG = newVectorWithInit(localMesh->ncoord);
    insertDirichlet(localSolCG, localMesh, u_D);
 
    // Save time for one process after setting up the problem data (rhs etc)
    MPI_Barrier(grid);
    if (rank ==0){
	TIME_SAVE(4,tv);
    }

    /* Create barrier, so that all process start solving at the same time.
       Thus, evaluate runtime of the solver independently */
    MPI_Barrier(grid);

    // Solve the problem using CG Solver and accumulate result
    nofIt = solve_cg(localMapping, sm_local, rhs, localSolCG, grid, 10e-10,
		      200);
    
    /* Set another barrier to wait for all processes to finish */
    MPI_Barrier(grid);
    // Save the total timevalue after all processes finished solving
    if (rank==0){
	TIME_SAVE(5, tv);
    }

    glblSolCG = accumulateResult(localMapping, localSolCG, grid);   
   
    // Save time after result is present
    MPI_Barrier(grid);
    if (rank==0){
	TIME_SAVE(6, tv);
    }
    

    if (rank ==0){
	delete2DMeshMapping(mapping, gridDims[0]);
	mesh_free(meshRefined);
	mesh_free(m1);
	free(localSolCG);

    } else{
	if (localMapping != NULL){
	    mesh_free(localMapping->localMesh);
	    deleteMeshMapping(localMapping);
	}
	free(glblSolCG);
    }



    return nofIt;
}


index getNumDegreesOfFreedom(int nrefines){
    index nNodes=2;
    index nDof=0;
    for (int i=2; i<=nrefines; ++i){
	nNodes*=2;
    }
    nNodes+=1;
    nDof=nNodes*nNodes;

    return nDof;

}

int main(int argc, char *argv[]){
    // Init MPI
    MPI_Init(&argc, &argv);
    
    char fname_p1[32] = "../Problem/problem1";
    char fname_result[40] =  "bench_results/bench_results_cg_sep.csv";
    FILE *fileResult; 
    struct timeval tv[7];
    index nDOF;
    int numRefines, nofIterations;

    // check input parameters
    if (argc < 4){
	fprintf(stderr, "Not enough input parameters!\n");
	abort();
    }

    // get grid dimensions and max refines from commandline arguments
    int gridDims[2] = {atoi(argv[1]), atoi(argv[2])};
    numRefines = atoi(argv[3]);


    // Initialize MPI grid Communication
    int nof_processes; MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int periods[2] = {false,false};

    MPI_Dims_create(nof_processes, 2, gridDims);
    MPI_Comm grid;
    MPI_Cart_create(MPI_COMM_WORLD, 2, gridDims, periods, true, &grid);
    MPI_Comm_rank(grid, &rank);

    // initialize result file
    if (rank == 0){
	printf("Init result file\n");
	initResultFile(fname_result, fileResult);
	fileResult = fopen(fname_result, "a");
    }

    if (rank==0){
	printf("=== Start bench_cg ===\n");
	printf("Argumens dimX=%d, dimY=%d, nRefines=%d\n", gridDims[0],gridDims[1],
		numRefines);
    }

    // solve poisson problem
    nofIterations = benchPoissonCG(fname_p1, numRefines,  F_vol, g_Neu, tv, grid);

    if (rank==0){
	nDOF = getNumDegreesOfFreedom(numRefines);
	printf("nDOF=%td\n",nDOF);
	printf("nofIterations=%d\n",nofIterations);
	saveBenchResults(tv, 7, nDOF, nofIterations, nof_processes, fileResult); 
	fclose(fileResult);
    }
    MPI_Finalize();

}
