#include <mpi.h>

#include "hpc.h"

/**
 * @brief Get an array with nbor ranks in the following order
 * 7 0 1
 * 6 x 2
 * 5 4 3
 * where x is the calling grid-node
 * @return int* Pointer to array of length 8 with nbor ranks
 */
int *getNbors(MPI_Comm grid) {
    int dims[2];
    int periods[2];
    int coords[2];
    MPI_Cart_get(grid, 2, dims, periods, coords);

    int *nbors = (int *)calloc(8, sizeof(int));
    if (!nbors) {
        printf("[E] Could not allocate memory in getNbors!\n");
        abort();
    }

    // dx and dy encode the x- and y-direction of the respective nbor
    int dx[8] = {0, 1, 1, 1, 0, -1, -1, -1};
    int dy[8] = {1, 1, 0, -1, -1, -1, 0, 1};
    for (int i = 0; i < 8; i++) {
        nbors[i] = MPI_PROC_NULL;
        int nborCoords[2] = {coords[0] + dx[i], coords[1] + dy[1]};
        if (nborCoords[0] >= 0 && nborCoords[0] < dims[0] &&
            nborCoords[1] >= 0 && nborCoords[1] < dims[1]) {
            MPI_Cart_rank(grid, nborCoords, &nbors[i]);
        }
    }
    return nbors;
}

/**
 * @brief
 *
 * @param vec
 * @param nbors
 * @param sendNborIndices
 * @param recvNborIndices
 * @param grid
 */
void accumulateVectorVComm(double *vec, int nbors[8],                       //
                           int sendNborIndices[4], int recvNborIndices[4],  //
                           MPI_Comm grid) {
    // Set up comm variables
    MPI_Request requests[8];
    int requestIndex = 0;
    double buffer[4];

    int rank;
    MPI_Comm_rank(grid, &rank);

    // Send
    for (int i = 0; i < 4; i++) {
        if (i == 0) {
            int coords[2];
            MPI_Cart_coords(grid, rank, 2, coords);
            printf("%d=(%d,%d) -> %d\n", rank, coords[0], coords[1],
                   nbors[sendNborIndices[i]]);
        }
        MPI_Isend(&vec[i], 1, MPI_DOUBLE,        //
                  nbors[sendNborIndices[i]], 0,  //
                  grid, &requests[requestIndex++]);
    }

    // Recv
    for (int i = 0; i < 4; i++) {
        if (i == 0) printf("%d <- %d\n", rank, nbors[recvNborIndices[i]]);
        MPI_Irecv(&buffer[i], 1, MPI_DOUBLE,     //
                  nbors[recvNborIndices[i]], 0,  //
                  grid, &requests[requestIndex++]);
    }

    // Wait for communication
    for (int i = 0; i < 8; i++) {
        // if (rank == 0) printf("%d\n", i);
        MPI_Status status;
        MPI_Wait(&requests[i], &status);
        // Ignore status for now
    }

    // Accumulate
    for (int i = 0; i < 4; i++) {
        vec[i] += buffer[i];
    }
}

void accumulateVectorV(double *vec, MPI_Comm grid) {
    // Get relevant nbors in grid
    int *nbors = getNbors(grid);

    /**
     * Sending directions of V nodes
     * up right
     * .-.
     * | |
     * .-.
     * left down
     *
     * Fixing a single distributed V node, we send
     * ne -> nw -> sw -> se -> ne.
     * This is a "circular" communication.
     */
    int rank;
    MPI_Comm_rank(grid, &rank);
    if (rank == 0) printf("round 1\n");
    int sendNborIndices1[4] = {6, 4, 0, 2};
    int recvNborIndices1[4] = {4, 2, 6, 0};
    accumulateVectorVComm(vec, nbors, sendNborIndices1, recvNborIndices1, grid);
    if (rank == 0) printf("round 2\n");

    /**
     * Fixing a single distributed V node, we send
     * ne <-> sw and nw <-> se.
     * This is a "cross" communication.
     */
    int nborIndices2[4] = {5, 3, 7, 1};
    accumulateVectorVComm(vec, nbors, nborIndices2, nborIndices2, grid);
}

/**
 * @brief Gets a distributed vector and changes it to an accumulated vector
 *
 * @param[in] mapping Mesh mapping
 * @param[in,out] vec Vector of size mapping->localMesh->ncoord. When calling
 * the function, this is a distributed vector. After the function terminates,
 * this is an accumulated vector.
 * @param[in] comm Comm for the accumulation
 *
 * @author Joselinda
 */
void accumulateVector(MeshMapping *mapping, double *vec, MPI_Comm comm) {
    mesh *localMesh = mapping->localMesh;
    index *vertexL2G = mapping->vertexL2G;
    index gncoord = mapping->globalNcoord;
    index lncoord = localMesh->ncoord;

    // Allocate globally sized, zero initialized vector
    // Add local vector entries to the corresponding global entries
    double *tmp = newVectorWithInit(gncoord);
    for (index i = 0; i < lncoord; i++) {
        tmp[vertexL2G[i]] += vec[i];
    }

    // Accumulate global vector
    MPI_Allreduce(MPI_IN_PLACE, tmp, gncoord, MPI_DOUBLE, MPI_SUM, comm);

    // Save local accumulated vector
    for (index i = 0; i < lncoord; i++) {
        vec[i] = tmp[vertexL2G[i]];
    }

    free(tmp);
}