#include <mpi.h>

#include "hpc.h"
#include <blas_level1.h>

/**
 * @brief Solve a system of linear equations distributed to multiple nodes.
 *
 * @param sm
 * @param x
 * @param y
 */
/*
void solve_cg(MeshMapping *localMapping, sed *localSM, double *rhs, double *u_local,
	      double *u_glbl, MPI_Comm grid, double tol, index maxIt)*/
void hpc_jacobi(MeshMapping *mapping, const sed *localSM, double *f,
                double *u, MPI_Comm grid, double tol, index maxIt) {

    double sigma, sigma_0;
    int rank, k;

    // get rank of node
    MPI_Comm_rank(grid, &rank);

    // accumulated and distributed vectors
    double *ww = newVectorWithInit(localSM->n);
    double *uu = newVectorWithInit(localSM->n);

    double *u0 = newVectorWithInit(localSM->n);

    // initialize vector for local inverse diagonal
    double *invDiag = newVectorWithInit(localSM->n);
    double *dd = newVectorWithInit(localSM->n);

    // initialze vector for right hand side and residual
    double *ff = newVectorWithInit(localSM->n);
    double *r = newVectorWithInit(localSM->n);

    // multiplied accumulated dd and ww -> dw
    double *dw = newVectorWithInit(localSM->n);

    // initially copy and accumulate diagonal of stiffness matrix
    blasl1_dcopy(localSM->x, dd, localSM->n, 1.0);
    accumulateVector(mapping, dd, grid);

    // invert diagonal of stiffness matrix and save for interation
    for (size_t i = 0; i < localSM->n; ++i)
        invDiag[i] = 1 / localSM->x[i];

    // set initial local u to the rhrs divided by the diagonal
    for (size_t i = 0; i < localSM->n; ++i)
        u0[i] = f[i] * invDiag[i];

    // sed_spmv_sym(localSM, u0, ff, -1.0, 1.0);

    // calculate residual and accumulate the vector
    sed_spmv_sym(localSM, u0, ff, 1.0, 0);
    blasl1_daxpy(ff, r, localSM->n, -1.0, 1.0);
    blasl1_dcopy(r, ww, localSM->n, 1.0);
    accumulateVector(mapping, ww, grid);

    sigma = dot_dist(ww, r, localSM->n, grid);
    sigma_0 = sigma;

    k = 0;
    // solve iteratively
    while (pow(tol, 2) * sigma_0 > sigma) {
        k += 1;

        // calculate local u
        for (size_t i = 0; i < localSM->n; ++i)
            dw[i] = dd[i] * ww[i];

        blasl1_daxpy(u0, uu, localSM->n, 1.0, 1.0);

        // calculate residual and accumulate the vector
        sed_spmv_sym(localSM, uu, ff, 1.0, 0);
        blasl1_daxpy(ff, r, localSM->n, -1.0, 1.0);
        blasl1_dcopy(r, ww, localSM->n, 1.0);
        accumulateVector(mapping, ww, grid);

        // calculate error
        sigma = dot_dist(ww, r, localSM->n, grid);
    }
}
