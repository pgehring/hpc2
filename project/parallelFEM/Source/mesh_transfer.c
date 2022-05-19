#include "hpc.h"


void transfer_MeshMetadata(index *buffer, MeshMapping ***mapping, int rankRecv,
			   int indXRecv, int indYRecv, MPI_Comm grid){
    int rank; MPI_Comm_rank(grid, &rank);

    int result = 1;
    
    if (rank==0){
	index sendData[5] = {mapping[indXRecv][indYRecv]->localMesh->ncoord,
			     mapping[indXRecv][indYRecv]->localMesh->nelem,
			     mapping[indXRecv][indYRecv]->localMesh->nedges,
			     mapping[indXRecv][indYRecv]->localMesh->nbdry,
			     mapping[indXRecv][indYRecv]->localMesh->nfixed};

	MPI_Send(sendData, 5, MPI_LONG_LONG, rankRecv, 1, grid);

    } else{
	MPI_Status *status;
	MPI_Recv(buffer, 5, MPI_LONG_LONG, 0, 1, grid, status);

    }
	
}

void transfer_MeshData(MeshMapping ***mapping, int rankRecv, mesh *localMesh,
		      int indXRecv, int indYRecv, MPI_Comm grid){
    
    int rank; MPI_Comm_rank(grid, &rank);

    if (rank=0){
	// get number of elements in each vector for sending
	index nCoord = mapping[indXRecv][indYRecv]->localMesh->ncoord;
	//index nElem = mapping[indXRecv][indYRecv]->localMesh->nelem;	 
	//index nBdry = mapping[indXRecv][indYRecv]->localMesh->nbdry;
	//index nFixed = mapping[indXRecv][indYRecv]->localMesh->nfixed;
	//index nEdges = mapping[indXRecv][indYRecv]->localMesh->nedges;
	//index nEdge2no = 2*nEdges;

	// Transfer coordinates
	printf("Sending meshdata to rank %d\n",rankRecv);
	MPI_Send(mapping[indXRecv][indYRecv]->localMesh->coord, nCoord,
		MPI_DOUBLE, rankRecv, 2, grid);
    } else{
	// Receive coordinates
	MPI_Status status;
	index recvCount;
	printf("Rank %d, retrieving mesh data\n",rankRecv);
	MPI_Recv(localMesh->coord, localMesh->ncoord, MPI_DOUBLE, 0, 2, grid,
		&status);

    }


}

mesh *mesh_transfer(MeshMapping ***mapping, MPI_Comm grid){
    int rank; MPI_Comm_rank(grid, &rank);

    // Get total number of processes
    int nof_processes; MPI_Comm_size(grid, &nof_processes);


    mesh *localMesh;
    index metadata[5];
    int proc_coords[2];

    if (rank == 0){
	// Do the transfer for all ranks
	for (int recvRank=1; recvRank<nof_processes; ++recvRank){
	    // Get grid coords for current rank
	    MPI_Cart_coords(grid, recvRank, 2, proc_coords);

	    printf("Process with rank %d has coords (%d,%d)\n",recvRank,
		   proc_coords[0], proc_coords[1]);
	
	    // Send metadata of the mesh from root to other processe
	    transfer_MeshMetadata(metadata, mapping, recvRank, proc_coords[0],
				  proc_coords[1], grid);

	    // get pointer to local mesh directly
	    localMesh = mapping[0][0]->localMesh;


	    printf("Going to transfer mesh data to Rank %d\n",recvRank);

	    // Send the mesh data
	    //transfer_MeshData(mapping, recvRank, localMesh, proc_coords[0],
	    //		      proc_coords[1], grid); 
	
	}

    } else{
	// Call the function with rank!=0 to get the metadata params 4 and 5
	// are not needed if rank!=0
	// metadata=[ncoord,nelem,nedges,nbdry,nfixed]
	transfer_MeshMetadata(metadata, mapping, rank, -1, -1, grid);

	// Allocate memory for local mesh
	localMesh = mesh_alloc(metadata[0],metadata[1],metadata[3]);

	// Receive the mesh data
	// params 4 and 5 are again not needed
	//transfer_MeshData(mapping, rank, localMesh, -1, -1, grid);
    }

    return localMesh;
}
