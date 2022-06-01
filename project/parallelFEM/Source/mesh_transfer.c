#include "hpc.h"

/**
  * @brief function to transfer the metadata, i.e the counts of vertices, elements
           etc. from root to the respective nonroot process. The metadata is needed
           to allocate memory for the mesh and mapping data from the nonroot process.
  * @param  *buffer buffer for the metadata on the receiving nonroot process.
            Can be nullpointer for the root process calling the function.
  * @param  ***mapping pointer to two dimensional array of MeshMapping structs i.e
            the global mapping stored on root. Can be nullpointer for nonroot
            processes calling the function
  * @param  rankRecv rank of the receiving nonroot process in the MPI grid
  * @param  indXRecv X-Index of the receiving nonroot process in the MPI grid
  * @param  indYRecv Y-Index of the receiving nonroot process in the MPI grid
  * @param  grid MPI communication domain organized in a two dimensional grid
*/
void transfer_Metadata(index *buffer, MeshMapping ***mapping, int rankRecv,
                       int indXRecv, int indYRecv, MPI_Comm grid) {
    int rank;
    MPI_Comm_rank(grid, &rank);

    int result = 1;

    if (rank == 0) {
        index sendData[8] = {mapping[indXRecv][indYRecv]->localMesh->ncoord,
                             mapping[indXRecv][indYRecv]->localMesh->nelem,
                             mapping[indXRecv][indYRecv]->localMesh->nedges,
                             mapping[indXRecv][indYRecv]->localMesh->nbdry,
                             mapping[indXRecv][indYRecv]->localMesh->nfixed,
                             mapping[0][0]->globalNcoord,
                             mapping[indXRecv][indYRecv]->lMeshDimX,
                             mapping[indXRecv][indYRecv]->lMeshDimY};

        DEBUG_PRINT("Sending metadata: [%td,%td,%td,%td,%td,%td,%td,%td]\n",
                    sendData[0], sendData[1], sendData[2], sendData[3],
                    sendData[4], sendData[5], sendData[6], sendData[7]);

        MPI_Send(sendData, 8, MPI_LONG_LONG, rankRecv, 1, grid);

    } else {
        MPI_Status status;
        MPI_Recv(buffer, 8, MPI_LONG_LONG, 0, 1, grid, &status);
    }
}

/**
  * @brief function to transfer the mesh data from root to the respective nonroot
           process. Since the local MappingStructs only contain pointers associated
           to the memory regions of the data, the data of each region have to be
           transferred seperately.
  * @param ***mapping pointer to two dimensional array of MeshMapping objects i.e
           the global mapping. Can be nullpointer for nonroot processes calling
           the function.
  * @param rankRecv rank of the receiving nonroot process in the MPI grid
  * @param *localMesh pointer to local MeshMapping struct of the nonroot process.
           Can be nullpointer for root process calling the function.
  * @param indXRecv, X-Index of the receiving nonroot process in the MPI grid
  * @param indYRecv, Y-Index of the receiving nonroot process in the MPI grid
  * @param grid MPI communication domain organized in a two dimensional grid
*/
void transfer_MeshData(MeshMapping ***mapping, int rankRecv, mesh *localMesh,
                       int indXRecv, int indYRecv, MPI_Comm grid) {
    int rank;
    MPI_Comm_rank(grid, &rank);

    if (rank == 0) {
        // get number of elements in each vector for sending
        index nCoord = mapping[indXRecv][indYRecv]->localMesh->ncoord;
        index nElem = mapping[indXRecv][indYRecv]->localMesh->nelem;
        index nBdry = mapping[indXRecv][indYRecv]->localMesh->nbdry;
        index nFixed = mapping[indXRecv][indYRecv]->localMesh->nfixed;
        index nEdge2no = mapping[indXRecv][indYRecv]->localMesh->nedges * 2;

        // Transfer coordinates
        DEBUG_PRINT("Sending meshdata to rank %d\n", rankRecv);
        // Send coordinates
        MPI_Send(mapping[indXRecv][indYRecv]->localMesh->coord, 2 * nCoord,
                 MPI_DOUBLE, rankRecv, 20, grid);

        // Send elements
        MPI_Send(mapping[indXRecv][indYRecv]->localMesh->elem, 7 * nElem,
                 MPI_LONG_LONG, rankRecv, 21, grid);

        // Send boundary
        MPI_Send(mapping[indXRecv][indYRecv]->localMesh->bdry, 4 * nBdry,
                 MPI_LONG_LONG, rankRecv, 22, grid);

        DEBUG_PRINT("Mesh data sent to rank %d\n", rankRecv);

    } else {
        // Receive coordinates
        MPI_Status status;
        index recvCount;
        DEBUG_PRINT("Rank %d, retrieving mesh data\n", rank);

        // Receive coordinates
        MPI_Recv(localMesh->coord, 2 * localMesh->ncoord, MPI_DOUBLE, 0, 20,
                 grid, &status);

        // Receive elements
        MPI_Recv(localMesh->elem, 7 * localMesh->nelem, MPI_LONG_LONG, 0, 21,
                 grid, &status);

        // Receive boundary
        MPI_Recv(localMesh->bdry, 4 * localMesh->nbdry, MPI_LONG_LONG, 0, 22,
                 grid, &status);

        DEBUG_PRINT("Rank %d, received mesh data!\n", rank);
    }
}

/**
  * @brief  function to transfer the mapping data, that is associated with the
            pointers in the local MeshMapping struct, from root to the respective
            process.
  * @param  ***globalMapping pointer to two dimensional array of MeshMapping structs
            i.e. the global mapping. Can be a nullpointer for nonroot processes
            calling the function
  * @param  *localMapping pointer to MappingStruct of the respective nonroot process.
            the associated memory region contains pointers to the local mapping data
            for the process. Can be nullpointer for the root process calling the function
  * @param  rankRecv rank of the process, the data shall be transferred to
  * @param  indXRecv X-index of the target process in the MPI grid
  * @param  indYRecv Y-index of the target process in the MPI grid
  * @param  grid MPI communication domain, i.e. the MPI grid
*/
void transfer_MappingData(MeshMapping ***globalMapping, MeshMapping *localMapping,
                          int rankRecv, int indXRecv, int indYRecv, MPI_Comm grid) {
    int rank;
    MPI_Comm_rank(grid, &rank);
    index ncoord, nelem, nbdry;

    // Send the data from root to process with rank=rankRecv
    if (rank == 0) {
        ncoord = globalMapping[indXRecv][indYRecv]->localMesh->ncoord;
        nelem = globalMapping[indXRecv][indYRecv]->localMesh->nelem;
        nbdry = globalMapping[indXRecv][indYRecv]->localMesh->nbdry;

        DEBUG_PRINT("Sending Mapping data to %d\n", rankRecv);
        // Transfer vertex mapping
        MPI_Send(globalMapping[indXRecv][indYRecv]->vertexL2G, ncoord,
                 MPI_LONG_LONG, rankRecv, 30, grid);

        // Transfer element mapping
        MPI_Send(globalMapping[indXRecv][indYRecv]->elemL2G, nelem,
                 MPI_LONG_LONG, rankRecv, 31, grid);

        // Transfer boundary mapping
        MPI_Send(globalMapping[indXRecv][indYRecv]->bdryL2G, nbdry,
                 MPI_LONG_LONG, rankRecv, 31, grid);

        DEBUG_PRINT("mapping sent to %d\n", rankRecv);

        // Receive mapping data on process with rank=rankRecv
    } else {
        index ncoord = localMapping->localMesh->ncoord;
        index nelem = localMapping->localMesh->nelem;
        index nbdry = localMapping->localMesh->nbdry;

        DEBUG_PRINT("rank %d receiving mappin data\n", rank);
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
        DEBUG_PRINT("rank %d received mapping data from root!\n", rank);
    }
}

/**
  * @brief function to transfer the local mapping objects with local meshes,
           stored in memory on the root process, to the other worker processes.
           Function can be called by the root and nonroot processes likewise. In
           the latter case, the parameter ***globalMapping can be a Nullpointer.
  * @param ***globalMapping pointer to two dimensional array of structs MeshMapping
           can be a Nullpointer for nonroot processes
  * @param grid two dimensional MPI grid, created using MPI_Cart_create()
  * @return pointer to local MeshMapping struct in the memory of the respective
            process.
*/
MeshMapping *mesh_transfer(MeshMapping ***globalMapping, MPI_Comm grid) {
    int rank;
    MPI_Comm_rank(grid, &rank);

    // Get total number of processes
    int nof_processes;
    MPI_Comm_size(grid, &nof_processes);

    MeshMapping *localMapping;
    index metadata[8];
    int proc_coords[2];
    mesh *localMesh;

    if (rank == 0) {
        // Do the transfer for all ranks
        for (int recvRank = 1; recvRank < nof_processes; ++recvRank) {
            // Get grid coords for current rank
            MPI_Cart_coords(grid, recvRank, 2, proc_coords);

            DEBUG_PRINT("Process with rank %d has coords (%d,%d)\n", recvRank,
                        proc_coords[0], proc_coords[1]);

            // Send metadata of the mesh from root to other processe
            transfer_Metadata(metadata, globalMapping, recvRank, proc_coords[0],
                              proc_coords[1], grid);

            DEBUG_PRINT("Going to transfer mesh data to Rank %d\n", recvRank);

            // Send the mesh data
            transfer_MeshData(globalMapping, recvRank, localMesh, proc_coords[0],
                              proc_coords[1], grid);

            transfer_MappingData(globalMapping, NULL, recvRank, proc_coords[0],
                                 proc_coords[1], grid);
        }

        // Get own coordinates
        MPI_Cart_coords(grid, 0, 2, proc_coords);
        localMapping = globalMapping[proc_coords[0]][proc_coords[1]];

    } else {
        // Call the function with rank!=0 to get the metadata. Params 4 and 5
        // are not needed if rank!=0
        // metadata=[ncoord,nelem,nedges,nbdry,nfixed, globalNCoords]
        transfer_Metadata(metadata, globalMapping, rank, -1, -1, grid);

        // Allocate memory for local mesh and insert data
        localMesh = mesh_alloc(metadata[0], metadata[1], metadata[3]);
        localMesh->nedges = metadata[2];
        localMesh->nfixed = metadata[4];

        DEBUG_PRINT("rank %d received Data: [%td,%td,%td,%td,%td,%td,%td,%td]\n", rank,
                    localMesh->ncoord, localMesh->nelem, localMesh->nedges,
                    localMesh->nbdry, localMesh->nfixed, metadata[5], metadata[6], metadata[7]);

        // Receive the mesh data
        // params 4 and 5 are again not needed
        transfer_MeshData(globalMapping, rank, localMesh, -1, -1, grid);

        // allocate mapping object
        localMapping = newMeshMapping(localMesh, metadata[5], metadata[6], metadata[7]);

        transfer_MappingData(NULL, localMapping, -1, -1, rank, grid);
    }

    return localMapping;
}
