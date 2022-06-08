#include <mpi.h>
#include "hpc.h"

/**
  * @brief Function to put the global result together. Collects elements from
           local part of accumulated vectors
  * @param localMapping pointer to local mapping struct
  * @param localResult pointer to double arrray of local result of a process
  * @param grid the MPI grid
  * @return the accumulated global solution
*/
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


/**
  * @brief function to insert the dirichlet data into a local solution vector.
	   used to initialize the vector for the solver. Vector should be 
	   zero initialized
  * @param u_local vector of local solution/approximation 
  * @param localMesh the local mesh of the calling process
  * @param u_D function pointer to functio that delivers value for the dirichtlet
	   boundary condition
*/
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


/**
  * @brief function to insert the indizes of dirichlet nodes, i.e. the nodes with
	   solution values fixed by dirichlet boundary conditions into the mesh
	   data structure. It's kinda a replacement of the function mesh_getFixed
	   that is applied to the solvers.
  * @mapping pointer to the local mapping of the process
*/
void mesh_getFixedNodes(MeshMapping *mapping){
    index nbdry = mapping->localMesh->nbdry;
    index *bdry = mapping->localMesh->bdry;

    // Allocate memory space with sufficient size
    index *fixed = calloc(2*nbdry,sizeof(index));
    index nfixed=0;

    /** Create array with flgs, stating whether a vertex has already
	been recognized as a fixed node*/
    int *setFlag = calloc(2*nbdry,sizeof(int));

    // Get fixed nodes
    for (index i=0; i<nbdry; ++i){
	if (!bdry[4*i+3]){
	    // Insert indices of fixed nodes if not already set. Update flags
	    if (!setFlag[bdry[4*i]]){
		fixed[nfixed] = bdry[4*i];
		setFlag[bdry[4*i]] = 1;
		nfixed++;
	    } 
	    if (!setFlag[bdry[4*i+1]]){
		fixed[nfixed] = bdry[4*i+1];
		setFlag[bdry[4*i+1]] = 1;
		nfixed++;
	    }
	}
    }

    // Write pointer and count to mesh
    mapping->localMesh->nfixed = nfixed;
    mapping->localMesh->fixed = fixed;
    
    free(setFlag);
}

/**
  * @brief function to set the values of a vector at the indizes, specified in
	   in localMesh->fixed to zero. Is needed for the solvers
  * @param localMesh pointer to the local mesh
  * @param x pointer to double array. Values at the respective indizes get 
	   overwritten
*/
void blockFixedNodes(mesh *localMesh, double *x){
    index nfixed = localMesh->nfixed;
    index *fixed = localMesh->fixed;

    for (index i=0; i<nfixed; ++i){
	x[fixed[i]] = 0;
    }

}
