#include "hpc.h"
#include <mpi.h>
#include "blas_level1.h"



/**
  * @brief Parallel CG-Solver for the poisson problem. Writes solution to u_local.
	   Needs to be called by all processes
  * @param localMapping pointer to the local mapping struct of the process
  * @param localSM pointer to local stiffness matrix as constructed by sed_sm_build
  * @param pointer to vector of local right hand side.
  * @param u_local pointer to double array for the local solution. Needs to be 
	   initialized with dirichlet boundary values at dirichlet indizes and 
	   0 elsewhere.
  * @param grid The MPI grid
  * @param tol Tolerance for the norm of the residuum
  * @param maxIt maximum number of iterations to be performed
  * @reutrn number of iterations performed by solver
*/
int solve_cg(MeshMapping *localMapping, sed *localSM, double *rhs, double *u_local,
	      MPI_Comm grid, double tol, index maxIt){
   
    int rank; MPI_Comm_rank(grid, &rank);
    index localDim =localMapping->localMesh->ncoord;
    int nofIt;

    // Declare and allocate necessary local distributed arrays
    double *r, *v;
    // Declare and allocate necessary local accumulated arrays
    double *w, *s;

    r = newVectorWithInit(localDim); 
    v = newVectorWithInit(localDim); 
    w = newVectorWithInit(localDim); 
    s = newVectorWithInit(localDim);

    // Declare scalars
    double sigma, sigmaOld, sigma0, alpha, beta, sv;

    int nofProcesses; MPI_Comm_size(grid, &nofProcesses);
    int debugRank=-1;

    // Update right hand side to accomodate for dirichlet data
    sed_spmv_sym(localSM, u_local, rhs, -1.0, 1.0);
   


    /* Initial computation of residual and its scalar product*/
    // Calculate matrix vectorprodukt Ku
    sed_spmv_sym(localSM, u_local, r, 1.0, 0); 
   
    // calculate residuum r=b-Ku
    blasl1_daxpy(r, rhs, localDim, 1, -1);   
    blockFixedNodes(localMapping->localMesh, r);
    
    // copy residuum for accumulation
    blasl1_dcopy(r, w, localDim, 1);	     

    // accumulate residuum
    accumulateVector(localMapping, w, grid);	     
    blockFixedNodes(localMapping->localMesh, w);


    // Copy in vector of search direction
    blasl1_dcopy(w, s, localDim, 1);		     
    // calculate scalar product of residuum
    sigma = dot_dist(w,r,localDim,grid);	     
   
    sigma0 = sigma;
    sigmaOld = sigma;

    /** Compute approximations iteratively */
    for (index iter=0; iter<maxIt; ++iter){
	nofIt = iter;

	// Compute spmv v=K*s and block fixed values in v
	sed_spmv_sym(localSM, s, v, 1.0, 0);        
	blockFixedNodes(localMapping->localMesh, v);
	
	// Compute step size
	sv = dot_dist(s,v,localDim,grid);
	alpha = sigma / sv;

	// Compute new approximation
	blasl1_daxpy(u_local, s, localDim, alpha, 1);
	
	// update residuum
	blasl1_daxpy(r, v, localDim, -alpha, 1);     

	// accumulate residuum
	blasl1_dcopy(r, w, localDim, 1);	    
	accumulateVector(localMapping, w, grid);
	blockFixedNodes(localMapping->localMesh, w);

	// calculate scalar product of residuum
	sigma = dot_dist(w,r,localDim,grid);	     
	beta = sigma/sigmaOld;

	sigmaOld = sigma;

	// calculate new search direction
	blasl1_daxpy(s, w, localDim, 1, beta);

	if (rank==1){
	    DEBUG_PRINT("sigma=%lf\n",sigma);
	}

	if (sqrt(sigma/sigma0)<=tol){
	    break;
	}
    }

    // Print output if naximum number of iterations has been reached
    if (nofIt == maxIt){
	fprintf(stderr, "hpc_solve_cg: Maximum number of iterations %td reached!\n",
		 maxIt);
    }

    // Free local arrays
    free(s);
    free(w);
    free(r);
    free(v);
    
    return nofIt;
}


