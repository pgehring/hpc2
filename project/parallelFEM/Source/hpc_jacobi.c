#include <mpi.h>

#include "hpc.h"
#include <blas_level1.h>

/**
 * @brief 
 * 
 * @param mapping 
 * @param localSM 
 * @param f 
 * @param u 
 * @param grid 
 * @param tol 
 * @param maxIt 
 */
int hpc_jacobi(MeshMapping *mapping, const sed *localSM, double *f,
                double *u, MPI_Comm grid, double tol, index maxIt) {

    double sigma, sigma_0;
    int rank, k;

    // get rank of node
    MPI_Comm_rank(grid, &rank);

    // accumulated and distributed vectors
    double *ww = newVectorWithInit(localSM->n);

    // initialize vector for local inverse diagonal
    double *invDiag = newVectorWithInit(localSM->n);
    double *dd = newVectorWithInit(localSM->n);

    // initialze vector for right hand side and residual
    double *r = newVectorWithInit(localSM->n);

    // multiplied accumulated dd and ww -> dw
    double *dw = newVectorWithInit(localSM->n);

    // update right hand side to accomodate for dirichlet data
    sed_spmv_sym(localSM, u, f, -1.0, 1.0);

    // initially copy and accumulate diagonal of stiffness matrix
    blasl1_dcopy(localSM->x, dd, localSM->n, 1.0);
    accumulateVector(mapping, dd, grid);

    // invert diagonal of stiffness matrix, save for interation
    // and apply to right hand side for initial value
    for (size_t i = 0; i < localSM->n; ++i){
        invDiag[i] = 1 / dd[i];
        u[i] = f[i] * invDiag[i];
    }

    // calculate residual and accumulate the vector
    sed_spmv_sym(localSM, u, r, 1.0, 0);
    blasl1_daxpy(r, f, localSM->n, 1.0, -1.0);
    blockFixedNodes(mapping->localMesh, r);

    blasl1_dcopy(r, ww, localSM->n, 1.0);

    // accumulate residuum
    accumulateVector(mapping, ww, grid);

    sigma = dot_dist(ww, r, localSM->n, grid);
    sigma_0 = sigma;

    k = 0;
    // solve iteratively
    while ((pow(tol, 2) * sigma_0 < sigma) && (k <= maxIt-1)) {
        k += 1;
        
        // calculate local u
        for (size_t i = 0; i < localSM->n; ++i)
            dw[i] = invDiag[i] * ww[i];

        blasl1_daxpy(u, dw, localSM->n, 1.0, 1.0);

        // calculate residual and accumulate the vector
        sed_spmv_sym(localSM, u, r, 1.0, 0);
        blasl1_daxpy(r, f, localSM->n, 1.0, -1.0);
        blockFixedNodes(mapping->localMesh, r);

        blasl1_dcopy(r, ww, localSM->n, 1.0);

        // accumulate residuum
        accumulateVector(mapping, ww, grid);

        // calculate error
        sigma = dot_dist(ww, r, localSM->n, grid);

    }
    return k;
}
