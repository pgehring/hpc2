#include "hpc.h"
#include <mpi.h>

double dot_dist(double *x_ac, double *y_dist, size_t len, MPI_Comm grid){
    int rank; MPI_Comm_rank(grid, &rank);
    int nof_processes; MPI_Comm_size(grid, &nof_processes);

    // Calculate local dot product
    double localDot = 0;
    double globalDot;
    for (size_t i=0; i<len; ++i){
	localDot += x_ac[i]*y_dist[i];
    }
    
    // Compute global result on root
    int res;
    res = MPI_Reduce(&localDot, &globalDot, 1, MPI_DOUBLE, MPI_SUM, 0, grid);
    if (res != MPI_SUCCESS){
	printf("Reduce operation for distributed dot product failed!\n");
	abort();
    }

    // Broadcast global result
    res = MPI_Bcast(&globalDot, 1, MPI_DOUBLE, 0, grid);
    if (res != MPI_SUCCESS){
	printf("Broadcast of distributed dot product failed!\n");
	abort();
    }

    return globalDot;
}


