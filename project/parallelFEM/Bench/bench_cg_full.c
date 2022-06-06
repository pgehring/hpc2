#include "hpc.h"
#include <mpi.h>
#include <time.h>
#include <sys/time.h>

/** Functions of neumann- and dirichlet boundary as well as the volume forces
  of Demo problem 1*/

struct timeval tv[50];
#define TIME_SAVE(j,tv) (gettimeofday(&tv[j], (struct timezone *)0))
#define TIME_ELAPSED(j, k, tv) \
	(1.E+3 * (tv[k].tv_sec - tv[j].tv_sec) + 1.E-3* (tv[k].tv_usec - tv[j].tv_usec))

double kappa(double x[2], index typ)
{
	return (1.0);
}

double u_D(double x[2])
{
	//  return ( 0.0 );
	return (x[0] * x[1]);
}

double randomEntry() { return sqrt(pow( ((double)rand()-RAND_MAX/2)/RAND_MAX,2.0));}

double F_vol(double x[2], index typ) { return (0.0); }

double g_Neu(double x[2], index typ) { return (x[0] * x[1]); }



void outputResults(struct timeval *tv, int nTimeStamps, index DOF, int nofIt, FILE *fileRes){
	double timesRel[nTimeStamps-1];	
	
	// Calculate times relative from start
	for (int i=1; i<nTimeStamps; ++i){
		//timesRel[i-1] = 1.E+3 *(double)(tv[i].tv_sec - tv[0].tv_sec) +
		//1.E-3*(double)(tv[i].tv_usec - tv[0].tv_usec);
		timesRel[i-1] = TIME_ELAPSED(0, i, tv);
	}
	
	// print results
	fprintf(fileRes,"%10td %10d ", DOF, nofIt);
	for (int i=0; i<nTimeStamps-2; ++i){
		fprintf(fileRes," %10lf", timesRel[i]);
	}
	fprintf(fileRes," %10lf\n",timesRel[nTimeStamps-2]);	

}


/** Wrapper function for solving the demo problem with test implementation
  of CG solver */
int solvePoissonCG(char *fname,  double (*fV)(double *, index), double (*fN)(double *, index),
		   int numRefines,  struct timeval *tv, int *dims, MPI_Comm grid){

	int rank; MPI_Comm_rank(grid, &rank);

	// Declare pointers for mapping structs
	MeshMapping ***mapping;
	mesh *m1;
	mesh *meshRefined;
	MeshMapping *localMapping;

	// Declare pointer to solution vectors and right hand side
	double *localSolCG, *glblSolCG, *rhs;

	// Number of needed iterations
	int maxit = 200;
	int nofIt;

	// Save initial Time
	TIME_SAVE(0, tv);

	// root: load and refine mesh, send local meshes
	if (rank==0){
		// Load basic mesh
		mesh *m1 = mesh_load(fname);
		mesh_getEdge2no(m1->nelem, m1->elem, &m1->nedges, &m1->edge2no);

		// initial refinement
		DEBUG_PRINT("Refine mesh\n","");
		meshRefined = mesh_initRefinement(m1, numRefines); 

		// create mapping
		DEBUG_PRINT("Create mapping\n","");
		mapping = mesh_split(meshRefined, dims);

		// transfer local meshes
		DEBUG_PRINT("Transfer meshes\n","");
		localMapping = mesh_transfer(mapping, grid);


		// Nonrot processes: Call method mesh_transfer to receive mesh and mapping
		// data
	} else{
		DEBUG_PRINT("Rank %d receivin local mesh\n",rank);
		localMapping = mesh_transfer(mapping, grid);
	}

	// extract fixed nodes of local mesh
	mesh_getFixedNodes(localMapping);

	// The local mesh
	mesh *localMesh = localMapping->localMesh;

	DEBUG_PRINT("Rank %d building local sm\n",rank);
	// Build the local stiffness matrices
	sed *sm_local = sed_sm_build(localMesh);

	// build local rhs
	DEBUG_PRINT("Rank %d building local rhs\n",rank);
	rhs = newVectorWithInit(localMesh->ncoord);
	mesh_build_rhs(localMesh, rhs, F_vol, g_Neu);

	// Synchronize processes to get timestamp where all have setup
	// the problem data T_Setup and save it. The solver can now also be evaluated
	// independently

	// Solve using CG Solver 
	// -------------------------------------------------
	// allocate and zero initialize result vector
	localSolCG = newVectorWithInit(localMesh->ncoord);
	// insert dirichlet data for indices of fixed nodes
	insertDirichlet(localSolCG, localMesh, u_D);

	DEBUG_PRINT("Rank %d calling CG solver\n",rank);
	// Solve the problem using CG Solver
	nofIt = solve_cg(localMapping, sm_local, rhs, localSolCG, grid, 10e-10,
			maxit);

	// Sychronize and save time T_Solve to get the timestamp where all processes
	// finished solving

	DEBUG_PRINT("Rank %d returned from solver, accumulating result\n",rank);
	// Accumulate the result (results only on root in actual global result)
	glblSolCG = accumulateResult(localMapping, localSolCG, grid);   

	// Save timestamp T_Res where the global result is present on root
	if (rank==0) TIME_SAVE(1, tv);
	

	//MPI_Barrier(grid);
	DEBUG_PRINT("Rank %d freeing memory\n",rank);
	// If rank=0 free the mapping, the refined mesh the base mesh and the
	if (rank ==0){
		delete2DMeshMapping(mapping, dims[0]);
		mesh_free(meshRefined);
		mesh_free(m1);
		free(localSolCG);
		free(glblSolCG);

		// if rank!=0 free the local mesh
	} else{
		mesh_free(localMapping->localMesh);
		deleteMeshMapping(localMapping);
		free(localSolCG);
	}

	DEBUG_PRINT("rank %d Returning from bench function\n",rank);
	return nofIt;
}

index getDegreesOfFreedom(int numRefines){
	index nodeCount = 2;
	index DOF=0;
	
	for (int i=1; i<numRefines; ++i){
		nodeCount *= 2;
	}
	nodeCount +=1;
	DOF = nodeCount*nodeCount;
	
	return DOF; 
} 

/** Bench the Workflow for solving with the CG Solver iteratively
 * 	  specified number of mesh refinements. Refinements, MPI grid dims
 * 	  and result file can be specified with the commandline arguments:
 * 	  1: xDim of MPI Grid
 * 	  2. yDim of MPI Grid
 * 	  3. number of mesh refinements to start with
 * 	  4. highest number of mesh refinements 
 * 	  5. filename to store the results in (Inside the directory bench_results)
*/

int main(int argc, char**argv){
	MPI_Init(&argc, &argv);

	// Initialize MPI grid Communication
	int nof_processes; MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int dims[2];
	int periods[2] = {false,false};

	// Maximum number of refines
	int numRefinesMax, numRefinesStart;	

	// File names of problem and results
	char fname_p1[32] = "../Problem/problem1";
	char *resDir = "bench_results/",  fname_res[64];	
	FILE *fileRes; 

	// Get grid dims ,number of refines, and result file from command line arguments
	if (argc < 6){
		fprintf(stderr, "Not enough input arguments\n");
		abort();
	} else {
		dims[0] = atoi(argv[1]);
		dims[1] = atoi(argv[2]);
		numRefinesStart = atoi(argv[3]);
		numRefinesMax = atoi(argv[4]);
		sprintf(fname_res, "%s%s", resDir, argv[5]);
	}
	
	// Create MPI grid
	MPI_Dims_create(nof_processes, 2, dims);
	MPI_Comm grid;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &grid);
	MPI_Comm_rank(grid, &rank);

	// Initialize timeval array for benchinig
	int numTimeVals = 2;
	struct timeval tv[numTimeVals];
	

	int nofIterations;
	index DOF;



	if (rank==0){
		printf("\n=== Start bench_cg_sep  ===\n","");
		printf("Grid dims: (%d,%d), refines start: %d, max refines: %d\n",
			dims[0],dims[1], numRefinesStart, numRefinesMax);
		printf("result file: %s\n",fname_res);	

		fileRes = fopen(fname_res,"a+");	
		fprintf(fileRes,"%10s %10s %10s \n", "nDOF","nIt", "T_Result"); 
	}

	MPI_Barrier(grid);

	for (int nrefine = numRefinesStart; nrefine <= numRefinesMax; ++nrefine){
		if (rank ==0){
			DEBUG_PRINT("\nNew iteration. Refines = %d\n",nrefine);
		}
		// Compute degrees of freedom
		DOF = getDegreesOfFreedom(nrefine);	
			
		DEBUG_PRINT("rank %d calling bench function\n",rank);
		
		// Call the function solving the poisson problem
		nofIterations = solvePoissonCG(fname_p1,  F_vol, g_Neu, nrefine, tv, dims,grid);
		
		DEBUG_PRINT("Rank %d finished solve function, nIter = %d\n",rank,nofIterations);
		
		// Output results and synchronize processes
		if (rank==0){
			outputResults(tv, numTimeVals, DOF, nofIterations, fileRes);
		}
		MPI_Barrier(grid);
	}

	if (rank==0){
		DEBUG_PRINT("\n== End bench_cg_sep ==\n","");
	}

	MPI_Finalize();
}











