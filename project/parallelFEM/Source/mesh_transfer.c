#include "hpc.h"


void transfer_Metadata(index *buffer, MeshMapping ***mapping, int rankRecv,
			   int indXRecv, int indYRecv, MPI_Comm grid){
    int rank; MPI_Comm_rank(grid, &rank);

    int result = 1;
    
    if (rank==0){
	index sendData[6] = {mapping[indXRecv][indYRecv]->localMesh->ncoord,
			     mapping[indXRecv][indYRecv]->localMesh->nelem,
			     mapping[indXRecv][indYRecv]->localMesh->nedges,
			     mapping[indXRecv][indYRecv]->localMesh->nbdry,
			     mapping[indXRecv][indYRecv]->localMesh->nfixed,
			     mapping[0][0]->globalNcoord};

	printf("Sending metadata: [%td,%td,%td,%td,%td,%td]\n",
		sendData[0], sendData[1], sendData[2], sendData[3],
		sendData[4], sendData[5]);

	MPI_Send(sendData, 6, MPI_LONG_LONG, rankRecv, 1, grid);

    } else{
	MPI_Status *status;
	MPI_Recv(buffer, 6, MPI_LONG_LONG, 0, 1, grid, status);

    }
	
}

void transfer_MeshData(MeshMapping ***mapping, int rankRecv, mesh *localMesh,
		      int indXRecv, int indYRecv, MPI_Comm grid){
    
    int rank; MPI_Comm_rank(grid, &rank);

    if (rank==0){
	// get number of elements in each vector for sending
	index nCoord = mapping[indXRecv][indYRecv]->localMesh->ncoord;
	index nElem = mapping[indXRecv][indYRecv]->localMesh->nelem;	 
	index nBdry = mapping[indXRecv][indYRecv]->localMesh->nbdry;
	index nFixed = mapping[indXRecv][indYRecv]->localMesh->nfixed;
	index nEdge2no = mapping[indXRecv][indYRecv]->localMesh->nedges*2;

	// Transfer coordinates
	printf("Sending meshdata to rank %d\n",rankRecv);
	// Send coordinates
	MPI_Send(mapping[indXRecv][indYRecv]->localMesh->coord, nCoord,
		MPI_DOUBLE, rankRecv, 20, grid);
	
	// Send elements
	MPI_Send(mapping[indXRecv][indYRecv]->localMesh->elem, nElem,
		MPI_LONG_LONG, rankRecv, 21, grid);

	// Send boundary
	MPI_Send(mapping[indXRecv][indYRecv]->localMesh->bdry, nBdry,
	MPI_LONG_LONG, rankRecv, 22, grid);

	printf("Mesh data sent to rank %d\n",rankRecv);

    } else{
	// Receive coordinates
	MPI_Status status;
	index recvCount;
	printf("Rank %d, retrieving mesh data\n",rank);
	
	// Receive coordinates 
	MPI_Recv(localMesh->coord, localMesh->ncoord, MPI_DOUBLE, 0, 20,
		grid,&status);
	
	// Receive elements
	MPI_Recv(localMesh->elem, localMesh->nelem, MPI_LONG_LONG, 0, 21,
	grid,&status);

	// Receive boundary
	MPI_Recv(localMesh->bdry, localMesh->nbdry, MPI_LONG_LONG, 0, 22,
	grid,&status);

	printf("Rank %d, received mesh data!\n",rank);

    }


}

void transfer_MappingData(MeshMapping ***globalMapping, MeshMapping *localMapping,
			  int rankRecv, int indXRecv, int indYRecv, MPI_Comm grid){

    int rank; MPI_Comm_rank(grid, &rank);
    index ncoord, nelem, nbdry;

    if (rank==0){
	ncoord = globalMapping[indXRecv][indYRecv]->localMesh->ncoord;
	nelem = globalMapping[indXRecv][indYRecv]->localMesh->nelem;
	nbdry = globalMapping[indXRecv][indYRecv]->localMesh->nbdry;

	printf("Sending Mapping data to %d\n",rankRecv);
	// Transfer verted mapping
	MPI_Send(globalMapping[indXRecv][indYRecv]->vertexL2G, ncoord,
		MPI_LONG_LONG, rankRecv, 30, grid);

	// Transfer element mapping
	MPI_Send(globalMapping[indXRecv][indYRecv]->elemL2G, nelem,
		MPI_LONG_LONG, rankRecv, 31, grid);

	// Transfer boundary mapping
	MPI_Send(globalMapping[indXRecv][indYRecv]->bdryL2G, nbdry,
		MPI_LONG_LONG, rankRecv, 31, grid);

	printf("mapping sent to %d\n",rankRecv);
    } else{
	index ncoord = localMapping->localMesh->ncoord;
	index nelem = localMapping->localMesh->nelem;
	index nbdry = localMapping->localMesh->nbdry;

	printf("rank %d receiving mappin data\n",rank);
	MPI_Status status;

	// Receive vertex mapping
	MPI_Recv(localMapping->vertexL2G, ncoord, MPI_LONG_LONG, 0, 30, grid,
		 &status);

	// Receive element mapping
	MPI_Recv(localMapping->elemL2G, nelem, MPI_LONG_LONG, 0, 31, grid,
		 &status);

	// Receive boundary mapping
	MPI_Recv(localMapping->bdryL2G, nbdry, MPI_LONG_LONG, 0, 31, grid,
		 &status);
	printf("rank %d received mapping data from root!\n",rank);
    }
}



MeshMapping *mesh_transfer(MeshMapping ***globalMapping, MPI_Comm grid){
    int rank; MPI_Comm_rank(grid, &rank);

    // Get total number of processes
    int nof_processes; MPI_Comm_size(grid, &nof_processes);


    MeshMapping *localMapping; 
    index metadata[6];
    int proc_coords[2];
    mesh *localMesh;
    
    if (rank == 0){
	// Do the transfer for all ranks
	for (int recvRank=1; recvRank<nof_processes; ++recvRank){
	    // Get grid coords for current rank
	    MPI_Cart_coords(grid, recvRank, 2, proc_coords);

	    printf("Process with rank %d has coords (%d,%d)\n",recvRank,
		   proc_coords[0], proc_coords[1]);
	
	    // Send metadata of the mesh from root to other processe
	    transfer_Metadata(metadata, globalMapping, recvRank, proc_coords[0],
	    		proc_coords[1], grid);

	    printf("Going to transfer mesh data to Rank %d\n",recvRank);

	    // Send the mesh data
	    transfer_MeshData(globalMapping, recvRank, localMesh, proc_coords[0],
	    		      proc_coords[1], grid);

	    transfer_MappingData(globalMapping, NULL, recvRank, proc_coords[0],
				 proc_coords[1], grid);
	
	}

	localMapping = globalMapping[0][0];

    } else{
	
	
	// Call the function with rank!=0 to get the metadata. Params 4 and 5
	// are not needed if rank!=0
	// metadata=[ncoord,nelem,nedges,nbdry,nfixed, globalNCoords]
	transfer_Metadata(metadata, globalMapping, rank, -1, -1, grid);

	// Allocate memory for local mesh and insert data
	localMesh = mesh_alloc(metadata[0],metadata[1],metadata[3]);
	localMesh->nedges = metadata[2];
	localMesh->nfixed = metadata[4];

	printf("rank %d received Data: [%td,%td,%td,%td,%td,%td]\n",rank,
	       localMesh->ncoord, localMesh->nelem, localMesh->nedges,
	       localMesh->nbdry,localMesh->nfixed,metadata[5] );

	// Receive the mesh data
	// params 4 and 5 are again not needed
	transfer_MeshData(globalMapping, rank, localMesh, -1, -1, grid);

	// allocate mapping object
	localMapping = newMeshMapping(localMesh, metadata[5]);

	transfer_MappingData(NULL, localMapping, -1, -1, rank, grid);
    }

    return localMapping;
}
