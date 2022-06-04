#include "hpc.h"
#include "blas_level1.h"
#include <time.h>

#ifndef MPI_GRID_X
#define MPI_GRID_X 2
#endif

#ifndef MPI_GRID_Y
#define MPI_GRID_Y 2
#endif

static size_t procUniqueShare = 5; // Number of vector elements that uniquely
				   // belong to one process

/**
  * @brief function to create random matrix entry
  * @return random double value
*/
double randomEntry() { return ((double)rand() - RAND_MAX / 2) / RAND_MAX *10; }

/**
  * @brief helper function to random initialize a vector
  * @param x double pointer to memory space for the vector
  * @param n dimension of the vector
*/
void vec_randomInit(double *x, index n){
    for (index i=0; i<n; ++i){
	x[i] = round(randomEntry()*100)/100;
    }
}

void vec_print(double *x, index n){
    for (index i=0; i<n; ++i){
	printf("%lf\n",x[i]);
    }
}
 
double *shareAccumulatedVector(double *vec, size_t n, MPI_Comm grid){
    int rank; MPI_Comm_rank(grid, &rank);
    int nof_processes; MPI_Comm_size(grid, &nof_processes);

    double *localVec = malloc((procUniqueShare+1)*sizeof(double));

    if (rank==0){
	size_t offset;
	// Share slices with overlapping elements
	for (int i=1; i<nof_processes-1; ++i){
	    DEBUG_PRINT("rank 0 sending slice of vector to rank %d\n",i);
	    offset = getSliceOffset(n, nof_processes, i);
	    MPI_Send(&vec[offset], procUniqueShare+1, MPI_DOUBLE,i, 1, grid);
	}
	
	// share slice non overlapping elements
	DEBUG_PRINT("rank 0 sending slice of vector to rank %d\n",nof_processes-1);
	offset = getSliceOffset(n, nof_processes, nof_processes-1);
	MPI_Send(&vec[offset], procUniqueShare, MPI_DOUBLE, nof_processes-1,
		 1, grid);

	// copy own slice
	DEBUG_PRINT("rank 0 copying own data\n");
	for (size_t i=0; i<procUniqueShare+1; ++i){
	    localVec[i] = vec[i];
	}
    } else if (rank != nof_processes-1){
	DEBUG_PRINT("rank %d receiving slice of vector from root\n",rank);
	MPI_Status status;
	MPI_Recv(localVec, procUniqueShare+1, MPI_DOUBLE, 0, 1, grid, &status);
    } else{
	DEBUG_PRINT("rank %d receiving slice of vector from root\n",nof_processes-1);
	MPI_Status status;
	MPI_Recv(localVec, procUniqueShare, MPI_DOUBLE, 0, 1, grid, &status);
    }

    return localVec;

}

double *shareDistributedVector(double *vec, size_t n, MPI_Comm grid){
    int rank; MPI_Comm_rank(grid, &rank);
    
    // share the vector
    double *localVec = shareAccumulatedVector(vec,  n, grid);

    if (rank==0){
	// account for accumulation of the shared vectors (addition of overlapping
	// elements)
	for (size_t i=1; i<n-1; ++i){
	    if (i%procUniqueShare == 0){
		vec[i] += vec[i]; 
	    }
	}
    }
    return localVec;
}


int main(int argc, char**argv){
    srand(time(NULL));
    
    MPI_Init(&argc, &argv);

    int dims[2] = {MPI_GRID_X, MPI_GRID_Y};
    int periods[2] = {false,false};
    int nof_processes; MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
    MPI_Dims_create(nof_processes, 2, dims);
    MPI_Comm grid; MPI_Cart_create(MPI_COMM_WORLD, 2 ,dims, periods, true, &grid);
    int rank; MPI_Comm_rank(grid, &rank);

    if (rank==0){
	printf("\n=== Start test_dot_dist ===\n");
    }

    size_t globalVecDim = nof_processes*procUniqueShare;
    size_t localVecDim;
    double *globalXVec, *globalYVec, *localXVec, *localYVec;
    double dotTest=0; double dotRef; double dotError=-1.0;

    // Determine dimensions of local vectors
    if (rank==nof_processes-1){
	localVecDim = procUniqueShare;
    } else{
	localVecDim = procUniqueShare+1;
    }

    // Allocate and initialize the global vectors
    if (rank==0){
	globalXVec = malloc(globalVecDim*sizeof(double));
	globalYVec = malloc(globalVecDim*sizeof(double));
	
	vec_randomInit(globalXVec, globalVecDim);
	vec_randomInit(globalYVec, globalVecDim);

	printf("\nReference y-vector:\n");
	vec_print(globalYVec, globalVecDim);
	printf("\nReference x-vector: \n");
	vec_print(globalXVec, globalVecDim);
    }

    // Receive local vectors
    localXVec = shareAccumulatedVector(globalXVec, globalVecDim, grid);
    localYVec = shareDistributedVector(globalYVec, globalVecDim, grid);
   
    // Compute the dot product with test implementation
    dotTest = dot_dist(localXVec, localYVec, localVecDim, grid); 

    // Compute the dot product with reference implementation and evaluate error
    if (rank==0){
	dotRef = blasl1_ddot(globalXVec,globalYVec,globalVecDim);
	dotError = dotRef - dotTest;
	printf("expected: %lf, actual: %lf, err=%lf\n",dotRef, dotTest,dotError);
    }

    // free all used memory
    free(localXVec);
    free(localYVec);
    
    if (rank==0){
        free(globalXVec);
        free(globalYVec);
        printf("\n=== End test_dot_dist ===\n");
    }

    MPI_Finalize();
}

