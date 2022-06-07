#include "hpc.h"
#include "blas_level1.h"
#include <mpi.h>
#include <sys/time.h>

#ifndef MPI_GRID_X
#define MPI_GRID_X 4
#endif
#ifndef MPI_GRID_Y
#define MPI_GRID_Y 2
#endif
#ifndef NUM_MESH_REFINES
#define NUM_MESH_REFINES 6
#endif

/** Functions of neumann- and dirichlet boundary as well as the volume forces
    of Demo problem 1*/

struct timeval tv[50];
#define TIME_SAVE(j) (gettimeofday(&tv[j], (struct timezone *)0))
#define TIME_ELAPSED(j, k) \
    (1.E+6 * (tv[k].tv_sec - tv[j].tv_sec) + (tv[k].tv_usec - tv[j].tv_usec))

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




/** Helper functions for debugging */
void printCoordinates(mesh *mesh){
    for (size_t i=0; i<mesh->ncoord; ++i){
	printf("%4zu: (%4lf,%4lf)\n", i, mesh->coord[i*2], mesh->coord[i*2+1]);
    }
}


void printMapping(MeshMapping *localMapping){
    for (index i=0; i<localMapping->localMesh->ncoord; ++i){
	printf("%td: %td\n",i, localMapping->vertexL2G[i]);
    }
}

void printBoundary(mesh *localMesh){
    for (index i=0; i<localMesh->nbdry; ++i){
	printf("bdry %td: (%td,%td), type:%td\n",
		localMesh->bdry[4*i+2], localMesh->bdry[4*i],
		localMesh->bdry[4*i+1], localMesh->bdry[4*i+3]);
    }
}

void index_print(index n, index *x){
    for (index i=0; i<n; ++i){
	printf("%td\n",x[i]);
    }
}

/** Wrapper function for solving the demo problem with test implementation
    of CG solver */
double *solvePoissonCG(char *fname, int numRefines, double (*fV)(double *, index),
		       double (*fN)(double *, index)){

    // Initialize MPI grid Communication
    int nof_processes; MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int dims[2] = {MPI_GRID_X,MPI_GRID_Y};
    int periods[2] = {false,false};

    MPI_Dims_create(nof_processes, 2, dims);
    MPI_Comm grid;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &grid);
    MPI_Comm_rank(grid, &rank);


    // Declare pointers for mapping structs
    MeshMapping ***mapping;
    mesh *m1;
    mesh *meshRefined;
    MeshMapping *localMapping;

    // Declare pointer to solution vectors and right hand side
    double *localSolCG, *glblSolCG, *rhs;
    index maxIt = 1.E+6;

    // root: load and refine mesh, send local meshes
    if (rank==0){
	
	// Load basic mesh
	mesh *m1 = mesh_load(fname);
	mesh_getEdge2no(m1->nelem, m1->elem, &m1->nedges, &m1->edge2no);
	
	// initial refinement
	DEBUG_PRINT("Refining mesh\n","");
	meshRefined = mesh_initRefinement(m1, numRefines); 
	DEBUG_PRINT("Count of nodes in new mesh: %td\n",meshRefined->ncoord);

	// create mapping
	mapping = mesh_split(meshRefined, dims);

	// transfer local meshes
	DEBUG_PRINT("Transferring mesh from root\n","");
	localMapping = mesh_transfer(mapping, grid);

    // Nonrot processes: Call method mesh_transfer to receive mesh and mapping
    // data
    } else{
	DEBUG_PRINT("rank %d Receiving mesh from root \n",rank);
	localMapping = mesh_transfer(mapping, grid);
    }
    // extract fixed nodes of local mesh
    mesh_getFixedNodes(localMapping);


    // The local mesh
    mesh *localMesh = localMapping->localMesh;

    // Build the local stiffness matrices
    DEBUG_PRINT("Rank %d building local stiffness matrix\n",rank);
    sed *sm_local = sed_sm_build(localMesh);

    // build local rhs
    rhs = newVectorWithInit(localMesh->ncoord);
    mesh_build_rhs(localMesh, rhs, F_vol, g_Neu);

    // Synchronize all processes ONLY FOR PRETTY COMMANDLINE OUTPUT
    int msg=0; MPI_Bcast(&msg, 1, MPI_INT, 0, grid);


    // Solve using CG Solver 
    // -------------------------------------------------
    // allocate and zero initialize result vector
    localSolCG = newVectorWithInit(localMesh->ncoord);
    // insert dirichlet data for indices of fixed nodes
    insertDirichlet(localSolCG, localMesh, u_D);
  
    // Solve the problem using CG Solver
    solve_cg(localMapping, sm_local, rhs, localSolCG, grid, 10e-10,
	     maxIt);

    // Accumulate the result (results only on root in actual global result)
    DEBUG_PRINT("rank %d: Accumulate Result\n", rank);
    glblSolCG = accumulateResult(localMapping, localSolCG, grid);   
   
    MPI_Barrier(grid);
    // If rank=0 free the mapping, the refined mesh the base mesh and the
    if (rank ==0){
	delete2DMeshMapping(mapping, dims[0]);
	mesh_free(meshRefined);
	//mesh_free(m1);
	free(localSolCG);

    // if rank!=0 free the local mesh
    } else{
	if (localMapping != NULL){
	    mesh_free(localMapping->localMesh);
	    deleteMeshMapping(localMapping);
	}
    }



    return glblSolCG;
}


// Wrapper function for solving the demo problem with test implementation
// of jacobi solver
double *solvePoissonJcb(char *fname, int numRefines, double (*fV)(double *, index),
		       double (*fN)(double *, index)){

    // Initialize MPI grid Communication
    int nof_processes; MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int dims[2] = {MPI_GRID_X,MPI_GRID_Y};
    int periods[2] = {false,false};

    MPI_Dims_create(nof_processes, 2, dims);
    MPI_Comm grid;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &grid);
    MPI_Comm_rank(grid, &rank);


    // Declare pointers for mapping structs
    MeshMapping ***mapping;
    mesh *m1;
    mesh *meshRefined;
    MeshMapping *localMapping;

    // Declare pointer to solution vectors and right hand side
    double *localSol, *glblSol, *rhs;
    index maxIt = 1.E+6;

    // root: load and refine mesh, send local meshes
    if (rank==0){
	
	// Load basic mesh
	mesh *m1 = mesh_load(fname);
	mesh_getEdge2no(m1->nelem, m1->elem, &m1->nedges, &m1->edge2no);
	
	// initial refinement
	DEBUG_PRINT("Refining mesh\n","");
	meshRefined = mesh_initRefinement(m1, numRefines); 
	DEBUG_PRINT("Count of nodes in new mesh: %td\n",meshRefined->ncoord);

	// create mapping
	mapping = mesh_split(meshRefined, dims);

	// transfer local meshes
	DEBUG_PRINT("Transferring mesh from root\n","");
	localMapping = mesh_transfer(mapping, grid);

    // Nonrot processes: Call method mesh_transfer to receive mesh and mapping
    // data
    } else{
	DEBUG_PRINT("rank %d Receiving mesh from root \n",rank);
	localMapping = mesh_transfer(mapping, grid);
    }

    // extract fixed nodes of local mesh
    mesh_getFixedNodes(localMapping);

    // The local mesh
    mesh *localMesh = localMapping->localMesh;

    // Build the local stiffness matrices
    DEBUG_PRINT("Rank %d building local stiffness matrix\n",rank);
    sed *sm_local = sed_sm_build(localMesh);

    DEBUG_PRINT("Rank %d building local rhs\n",rank);
    // build local rhs
    rhs = newVectorWithInit(localMesh->ncoord);
    mesh_build_rhs(localMesh, rhs, F_vol, g_Neu);

    // Synchronize all processes ONLY FOR PRETTY COMMANDLINE OUTPUT
    int msg=0; MPI_Bcast(&msg, 1, MPI_INT, 0, grid);


    // Solve using jacobi solver
    // -------------------------------------------------
  
    // allocate and zero initialize result vector
    localSol = newVectorWithInit(localMesh->ncoord);
    // insert dirichlet data for indices of fixed nodes
    insertDirichlet(localSol, localMesh, u_D);
  
    // Solve the problem using jacobi Solver
    DEBUG_PRINT("Rank %d Calling jacobi solver\n",rank);
    hpc_jacobi(localMapping, sm_local, rhs, localSol, grid, 10e-10, maxIt);

    // Accumulate the result (results only on root in actual global result)
    DEBUG_PRINT("rank %d: Accumulate Result\n",rank);
    glblSol = accumulateResult(localMapping, localSol, grid);   

    MPI_Barrier(grid);
    // If rank=0 free the mapping, the refined mesh the base mesh and the
    if (rank ==0){
	delete2DMeshMapping(mapping, dims[0]);
	mesh_free(meshRefined);
	//mesh_free(m1);
	free(localSol);

    // if rank!=0 free the local mesh
    } else{
	
	if (localMapping != NULL){
	    mesh_free(localMapping->localMesh);
	    deleteMeshMapping(localMapping);
	}
	
    }

    return glblSol;
}

/** Wrapper for solving the demo problem with the reference implementation */
double *solvePoissonRef(char *fname, int numRefines, double (*fV)(double *, index),
		     double (*fN)(double *, index)){

    index n, k, ncoord, nelem, nbdry, nfixed, nedges, total, *bdry,
        cnt = 0, N = 0, MAX_ITER = 100;
    double *b, *x, *w, *Coord, x1[2], x2[2], m[2];
    mesh **H, *T;
    sed **A;

    TIME_SAVE(0);
    N = numRefines;
    DEBUG_PRINT("Load data form %s, no. refinements = %g\n", fname, (double)N);

    /* Allocate memory for hierachy */
    H = (mesh **)malloc((N + 1) * sizeof(mesh *));
    A = (sed **)malloc((N + 1) * sizeof(sed *));

    /* Load problem */
    H[0] = mesh_load(fname); /* load geometry */
    mesh_getEdge2no(H[0]->nelem, H[0]->elem, &H[0]->nedges, &H[0]->edge2no);
    H[0]->fixed =
        mesh_getFixed(H[0]->ncoord, H[0]->bdry, H[0]->nbdry, &H[0]->nfixed);
    DEBUG_PRINT("\nInit mesh  # dofs =  %10g\n",
           (double)H[0]->ncoord + H[0]->nedges);


    /* Build stiffness matrix, refine mesh and create hierachy  */
    k = 0;
    while (1)
    {
        A[k] = sed_nz_pattern(H[k]); /* get pattern of matrix */
        if (!A[k])
            return NULL;
        if (!sed_buildS(H[k], A[k]))
            return NULL; /* assemble coefficient matrix */
        if (k >= N)
            break;
        H[k + 1] = mesh_refine(H[k]);
        mesh_getEdge2no(H[k + 1]->nelem, H[k + 1]->elem, &H[k + 1]->nedges,
                        &H[k + 1]->edge2no);
        H[k + 1]->fixed = mesh_getFixed(H[k + 1]->ncoord, H[k + 1]->bdry,
                                        H[k + 1]->nbdry, &H[k + 1]->nfixed);
        k++;
    }
    TIME_SAVE(1);
    printf("Final mesh # dofs =  %10g\n", (double)H[N]->ncoord+H[N]->nedges);
    printf("# refinements     =  %10g\n", (double)N);

    n = A[N]->n;
    x = (double *)calloc(n, sizeof(double)); /* get workspace for sol*/
    w = (double *)calloc(n, sizeof(double)); /* get temporary workspace */
    b = (double *)calloc(n, sizeof(double)); /* get workspace for rhs*/
    mesh_build_rhs(H[N], b, fV,
                  fN); /* build rhs (volume and Neumann data */
    TIME_SAVE(2);

    /* incorporate Dirichlet data */
    ncoord = H[N]->ncoord;
    nelem = H[N]->nelem;
    nbdry = H[N]->nbdry;
    bdry = H[N]->bdry;
    H[N]->fixed = mesh_getFixed(ncoord, bdry, nbdry, &H[N]->nfixed);
    nfixed = H[N]->nfixed;
    nedges = H[N]->nedges;
    Coord = H[N]->coord;

    for (k = 0; k < nbdry; k++)
    {
        if (!bdry[4 * k + 3])
        {
            x1[0] = Coord[2 * bdry[4 * k]];
            x1[1] = Coord[2 * bdry[4 * k] + 1];
            x2[0] = Coord[2 * bdry[4 * k + 1]];
            x2[1] = Coord[2 * bdry[4 * k + 1] + 1];
            m[0] = (x1[0] + x2[0]) / 2.0;
            m[1] = (x1[1] + x2[1]) / 2.0;
            x[bdry[4 * k]] = u_D(x1);
            x[bdry[4 * k + 1]] = u_D(x2);
            x[ncoord + bdry[4 * k + 2]] = u_D(m);
        }
    }

    TIME_SAVE(3);
    cnt = hpc_mg(A, b, x, 1e-10, 50, H, N, 2, 2, 1);
    TIME_SAVE(4);

    for (k = 0; k < HPC_MIN(10, A[N]->n); k++)
    {
        printf(" x[%g] = %g\n", (double)k, x[k]);
    }

    DEBUG_PRINT("\n","");
    DEBUG_PRINT("Time load & create hierarchy = %9i ns\n", (int)TIME_ELAPSED(0, 1));
    DEBUG_PRINT("Time building rhs            = %9i ns\n", (int)TIME_ELAPSED(1, 2));
    DEBUG_PRINT("Time Dirichlet values        = %9i ns\n", (int)TIME_ELAPSED(2, 3));
    DEBUG_PRINT("Time solve LSE               = %9i ns\n", (int)TIME_ELAPSED(3, 4));
    DEBUG_PRINT("No. iterations               = %9g\n", (double)cnt);
    DEBUG_PRINT("========================================\n\n","");

    DEBUG_PRINT("\nMemory\n","");
    DEBUG_PRINT("Coordinates : %12zu Byte\n", ncoord * 2 * sizeof(double));
    DEBUG_PRINT("Elements :    %12zu Byte\n", nelem * 7 * sizeof(index));
    DEBUG_PRINT("Boundary :    %12zu Byte\n", nbdry * 4 * sizeof(index));
    DEBUG_PRINT("Edge2no :     %12zu Byte\n", nedges * 2 * sizeof(index));
    total = ncoord * 2 * sizeof(double) +
            (7 * nelem + 4 * nbdry + nedges * 2) * sizeof(index);
    DEBUG_PRINT("Total :       %12.6g MByte\n", (double)total / 1024. / 1024.);

    for (k = 0; k <= N; k++)
    {
        mesh_free(H[k]);
        sed_free(A[k]);
    }
    free(H);
    free(A);

    return x;
}




/** Helper functions for the test evaluation */

index getSolutionDimension(int numMeshRefines){
    index n=2;
    for (int i=0; i<numMeshRefines; ++i){
	n += (index)pow(2.0,(double)i);
    }
    return n*n;
}

double computeRMSE(size_t n, double *x, double *y){
    double rmse=0;
    if (n == 0||!x ||!y){
	printf("NULL\n");
	return -1.0;
    }

    for (size_t i=0; i<n; ++i){
	rmse += pow(x[i]-y[i],2.0);
    }
    rmse = sqrt(rmse);
    return rmse;
}





int main(int argc, char**argv){
    MPI_Init(&argc, &argv);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nof_processes; MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);

    size_t dimSolution = getSolutionDimension(NUM_MESH_REFINES);

    if (rank==0){
	printf("\n=== Start test_solver ===\n");
	printf("-- number of degrees of freedom: %zu\n",dimSolution);
	printf("-- number of processes used: %d\n",nof_processes);
    }

    char fname_p1[32] = "../Problem/problem1";
   
    double *refSolP1;
    /* Solve using reference implementation. One refinement of the mesh less,
       sine a implicit refinement is performed during the creation of the 
       stiffness matrix. ->comparability to our solvers */
    if (rank==0){
	printf("\n- Solve problem using reference implementation...\n");
	refSolP1 = solvePoissonRef(fname_p1, NUM_MESH_REFINES-1, F_vol, g_Neu);
	printf("\nReference solution [0-5]: \n");
	vecPrint(refSolP1, 5);		
    }

    /* solve using jacobi and CG solver */
    double *cgSolP1, cgRMSEP1, *jcbSolP1, jcbRMSEP1;
    
    if (rank==0){
	printf("\n- Solve using CG solver...\n");
    }
    DEBUG_PRINT("rank %d calling CG solver\n",rank);
    
    cgSolP1 =  solvePoissonCG(fname_p1, NUM_MESH_REFINES, F_vol, g_Neu);
    
    DEBUG_PRINT("rank %d finished CG solver\n",rank);




    if (rank==0){
	printf("\n-  Solve using Jacobi solver...\n"); 
    }
    DEBUG_PRINT("rank %d calling Jacobi solver\n",rank);
    
    jcbSolP1 = solvePoissonJcb(fname_p1, NUM_MESH_REFINES, F_vol, g_Neu);
    
    DEBUG_PRINT("rank %d finished jacobi solver\n",rank);

    if (rank==0){
	printf("\nSolution using CG [0-5]:\n");
	vecPrint(cgSolP1, 5);
	cgRMSEP1 = computeRMSE(dimSolution, refSolP1, cgSolP1);
	printf("RMSE = %lf\n", cgRMSEP1);

	printf("\nSolution using w-Jacobi [0-5]:\n");
	vecPrint(jcbSolP1, 5);
	jcbRMSEP1 = computeRMSE(dimSolution, refSolP1, jcbSolP1);
	printf("RMSE = %lf\n", jcbRMSEP1);
    }


    DEBUG_PRINT("rank %d finalizing\n",rank);
    
    // free(refSolP1);
    // free(cgSolP1);


    MPI_Finalize();
}











