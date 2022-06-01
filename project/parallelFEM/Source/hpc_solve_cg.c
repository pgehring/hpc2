#include "hpc.h"
#include <mpi.h>
#include "blas_level1.h"

void vector_print(double *x,index len){
    if (!x || len==0){
	printf("NULL\n");
	return;
    }
    for (index i=0; i<len; ++i){
	printf("%td : %lf\n",i,x[i]);
    }
}

void solve_cg(MeshMapping *localMapping, sed *localSM, double *rhs, double *u_local,
	      double *u_glbl, MPI_Comm grid, double tol, index maxIt){
   
    int rank; MPI_Comm_rank(grid, &rank);
    index localDim =localMapping->localMesh->ncoord;

    // Declare and allocate necessary local distributed arrays
    double *r, *v;
    // Declare and allocate necessary local accumulated arrays
    double *w, *s;

    r = newVectorWithInit(localDim); 
    v = newVectorWithInit(localDim); 
    w = newVectorWithInit(localDim); 
    s = newVectorWithInit(localDim);

    // Declare scalars
    double sigma, sigmaOld, alpha, beta, sv;

    int nofProcesses; MPI_Comm_size(grid, &nofProcesses);


    int debugRank=-1;

    // Initial computation of residual and its scalar product
    
    // Calculate matrix vectorprodukt Ku
    sed_spmv_sym(localSM, u_local, r, 1.0, 0); 
   
    // calculate residuum
    blasl1_daxpy(r, rhs, localDim, 1, -1);   

    // copy residuum for accumulation
    blasl1_dcopy(r, w, localDim, 1);	     

    // accumulate residuum
    accumulateVector(localMapping, w, grid);	     
   
    // Copy in vector of search direction
    blasl1_dcopy(w, s, localDim, 1);		     
    // calculate scalar product of residuum
    sigma = dot_dist(w,r,localDim,grid);	     
    
    sigmaOld = sigma;

    /** Compute approximations iteratively */
    for (index iter=0; iter<maxIt; ++iter){
	// Compute spmv K*s
	sed_spmv_sym(localSM, s, v, 1.0, 0);        

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

	// calculate scalar product of residuum
	sigma = dot_dist(w,r,localDim,grid);	     
	beta = sigma/sigmaOld;

	sigmaOld = sigma;


	// calculate new search direction
	blasl1_daxpy(s, w, localDim, 1, beta);

	if (rank==1){
	    printf("sigma=%lf\n",sigma);
	}
    }

    // Free local arrays
    free(s);
    free(w);
    free(r);
    free(v);
}


