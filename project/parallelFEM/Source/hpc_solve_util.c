#include <mpi.h>
#include "hpc.h"


double *accumulateResult(MeshMapping *localMapping, double *localResult,
			 MPI_Comm grid){
    double *glblResult;
    index globalDim = localMapping->globalNcoord;
    index localDim = localMapping->localMesh->ncoord;

    // Allocate zero initialized global sized buffer
    double *buffer = newVectorWithInit(globalDim);
    
    // Get MPI params
    int rank; MPI_Comm_rank(grid, &rank);
    int nof_processes; MPI_Comm_size(grid, &nof_processes);

    // accumulate Result on root
    if (rank==0){
	// Allocate and zero initialize global solution
	glblResult = newVectorWithInit(globalDim);

	// Iterate over all nonroot processes
	for (int r=1; r<nof_processes; ++r){
	    // Receive local solution of current process 
	    MPI_Status status;
	    
	    DEBUG_PRINT("Receiving result from rank %d\n",r);
	    MPI_Recv(buffer, globalDim, MPI_DOUBLE, r, r+1, grid, &status);
	    DEBUG_PRINT("Result from %d received!\n",r);

	    // Write nonzero entries to global solution
	    for (index i=0; i<globalDim; ++i){
		if (buffer[i] != 0){
		    glblResult[i] = buffer[i];
		}
	    }
	}

	// Write entries of own local solution
	for (index i=0; i<localDim; ++i){
	    if (localResult[i] != 0){
		glblResult[localMapping->vertexL2G[i]] = localResult[i];
	    }
	}
    // nonroot processes: send elements of local solution in global context
    } else{
	// Write elements of local solution to right indizes in buffer
	for (index i=0; i<localDim; ++i){
	    buffer[localMapping->vertexL2G[i]] = localResult[i];
	}

	// Send the buffer to root
	DEBUG_PRINT("rank %d Sending result to root\n",rank);
	MPI_Send(buffer, globalDim, MPI_DOUBLE, 0, rank+1, grid);
    }

    free(buffer);

    return glblResult;
    
}


/** function to insert dirichlet boundary data in solution vector */
void insertDirichlet(double *u_local, mesh *localMesh, double (*u_D)(double *)){
    index nbdry = localMesh->nbdry;
    index nodeInds[2];
    double coords[2];

    // Search for dirichlet boundaries
    for (index i=0; i<nbdry; ++i){
	if (!localMesh->bdry[4*i+3]){
	    // Insert node indices of current boundary in vector
	    nodeInds[0] = localMesh->bdry[4*i];
	    nodeInds[1] = localMesh->bdry[4*i+1];

	    // Get the boundary value for each node and write it to the
	    // solution vector
	    for (int k=0; k<2; ++k){
		coords[0] = localMesh->coord[2*nodeInds[k]];	;
		coords[1] = localMesh->coord[2*nodeInds[k]+1];
		u_local[nodeInds[k]] = u_D(coords);
	    }
	}
    }

}
