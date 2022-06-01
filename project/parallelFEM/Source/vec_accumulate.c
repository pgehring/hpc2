#include <mpi.h>

#include "hpc.h"

#define NBOR_SW 0
#define NBOR_S 1
#define NBOR_SE 2
#define NBOR_W 3
#define NBOR_E 4
#define NBOR_NW 5
#define NBOR_N 6
#define NBOR_NE 7

/**
 * @brief Get an array with nbor ranks in the following order
 * 5 6 7
 * 3 x 4
 * 0 1 2
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
    int dx[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
    int dy[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    for (int i = 0; i < 8; i++) {
        nbors[i] = MPI_PROC_NULL;
        int nborCoords[2] = {coords[0] + dx[i], coords[1] + dy[i]};
        if (nborCoords[0] >= 0 && nborCoords[0] < dims[0] &&
            nborCoords[1] >= 0 && nborCoords[1] < dims[1]) {
            MPI_Cart_rank(grid, nborCoords, &nbors[i]);
        }
    }
    return nbors;
}

/**
 * @param grid
 * @return if the node is on the north border on the comm grid
 */
bool onNorthBorder(MPI_Comm grid) {
    int dims[2];
    int periods[2];
    int coords[2];
    MPI_Cart_get(grid, 2, dims, periods, coords);
    return coords[1] + 1 == dims[1];
}

/**
 * @param grid
 * @return if the node is on the east border on the comm grid
 */
bool onEastBorder(MPI_Comm grid) {
    int dims[2];
    int periods[2];
    int coords[2];
    MPI_Cart_get(grid, 2, dims, periods, coords);
    return coords[0] + 1 == dims[0];
}

/**
 * @param grid
 * @return if the node is on the south border on the comm grid
 */
bool onSouthBorder(MPI_Comm grid) {
    int dims[2];
    int periods[2];
    int coords[2];
    MPI_Cart_get(grid, 2, dims, periods, coords);
    return coords[1] == 0;
}

/**
 * @param grid
 * @return if the node is on the north border on the comm grid
 */
bool onWestBorder(MPI_Comm grid) {
    int dims[2];
    int periods[2];
    int coords[2];
    MPI_Cart_get(grid, 2, dims, periods, coords);
    return coords[0] == 0;
}

/**
 * @brief This function sends vec[0..3] to the nbors specified in
 * sendNborIndices, receives buf[0..3] values from the nbors specified
 * in recvNborIndices, and adds them.
 *
 * @param vec vec with V nodes as the first 4 entries. Must be ordered like this
 * in the mesh:
 * 2 3
 * 0 1
 * @param nbors Array with nbor ranks as returned by getNbors
 * @param sendNborIndices
 * @param recvNborIndices
 * @param grid Comm grid
 */
void accumulateVectorVComm(double *vec, int nbors[8],                       //
                           int sendNborIndices[4], int recvNborIndices[4],  //
                           MPI_Comm grid) {
    // Set up comm variables
    MPI_Request requests[8];
    int requestIndex = 0;
    double buffer[4] = {0, 0, 0, 0};

    // Send
    for (int i = 0; i < 4; i++) {
        MPI_Isend(&vec[i], 1, MPI_DOUBLE,        //
                  nbors[sendNborIndices[i]], 0,  //
                  grid, &requests[requestIndex++]);
    }

    // Recv
    for (int i = 0; i < 4; i++) {
        MPI_Irecv(&buffer[i], 1, MPI_DOUBLE,     //
                  nbors[recvNborIndices[i]], 0,  //
                  grid, &requests[requestIndex++]);
    }

    // Wait for communication
    for (int i = 0; i < 8; i++) {
        MPI_Status status;
        MPI_Wait(&requests[i], &status);
        // Ignore status for now
    }

    // Accumulate
    for (int i = 0; i < 4; i++) {
        vec[i] += buffer[i];
    }
}

/**
 * @brief Accumulates the V nodes.
 *
 * @param vec The V node entries must be the first four values and
 * they must be ordered like this in the mesh:
 * 2 3
 * 0 1
 * @param grid Comm grid
 */
void accumulateVectorV(double *vec, MPI_Comm grid, int nbors[8]) {
    /**
     * Fixing a single distributed V node, we send
     * n -> w -> s -> e -> n.
     * This is a "circular" communication.
     */
    int sendNborIndices1[4] = {onWestBorder(grid) && !onSouthBorder(grid) ? NBOR_S : NBOR_W,  // Nodes at the border
                               onSouthBorder(grid) && !onEastBorder(grid) ? NBOR_E : NBOR_S,  // can not communicate
                               onNorthBorder(grid) && !onWestBorder(grid) ? NBOR_W : NBOR_N,  // circularly :(
                               onEastBorder(grid) && !onNorthBorder(grid) ? NBOR_N : NBOR_E};
    int recvNborIndices1[4] = {onSouthBorder(grid) && !onWestBorder(grid) ? NBOR_W : NBOR_S,  //
                               onEastBorder(grid) && !onSouthBorder(grid) ? NBOR_S : NBOR_E,  //
                               onWestBorder(grid) && !onNorthBorder(grid) ? NBOR_N : NBOR_W,  //
                               onNorthBorder(grid) && !onEastBorder(grid) ? NBOR_E : NBOR_N};
    accumulateVectorVComm(vec, nbors, sendNborIndices1, recvNborIndices1, grid);

    /**
     * Fixing a single distributed V node, we send
     * ne <-> sw and nw <-> se.
     * This is a "cross" communication.
     * After that, all nodes have the same value for V.
     */
    int nborIndices2[4] = {NBOR_SW, NBOR_SE, NBOR_NW, NBOR_NE};
    accumulateVectorVComm(vec, nbors, nborIndices2, nborIndices2, grid);
}

void accumulateVectorE(MeshMapping *mapping, double *vec, MPI_Comm grid, int nbors[8]) {
    index lMeshDimX2 = mapping->lMeshDimX - 2;
    index lMeshDimY2 = mapping->lMeshDimY - 2;

    MPI_Request requests[8];
    int requestIndex = 0;
    index offsets[4] = {4,                   //
                        4 + lMeshDimX2,      //
                        4 + 2 * lMeshDimX2,  //
                        4 + 2 * lMeshDimX2 + lMeshDimY2};
    index counts[4] = {lMeshDimX2, lMeshDimX2, lMeshDimY2, lMeshDimY2};

    int sendNborIndices[4] = {NBOR_S, NBOR_N, NBOR_W, NBOR_E};
    for (int i = 0; i < 4; i++) {
        MPI_Isend(&vec[offsets[i]], counts[i], MPI_DOUBLE,  //
                  nbors[sendNborIndices[i]], 0,             //
                  grid, &requests[requestIndex++]);
    }

    index nofENodesP4 = 4 + 2 * lMeshDimX2 + 2 * lMeshDimY2;
    double buffer[nofENodesP4];
    for (int i = 0; i < nofENodesP4; i++) buffer[i] = 0;

    int recvNborIndices[4] = {NBOR_N, NBOR_S, NBOR_E, NBOR_W};
    offsets[0] = 4 + lMeshDimX2;
    offsets[1] = 4;
    offsets[2] = 4 + 2 * lMeshDimX2 + lMeshDimY2;
    offsets[3] = 4 + 2 * lMeshDimX2;
    for (int i = 0; i < 4; i++) {
        MPI_Irecv(&buffer[offsets[i]], counts[i], MPI_DOUBLE,  //
                  nbors[recvNborIndices[i]], 0,                //
                  grid, &requests[requestIndex++]);
    }

    for (int i = 0; i < 8; i++) {
        MPI_Status status;
        MPI_Wait(&requests[i], &status);
    }

    for (int i = 4; i < nofENodesP4; i++) {  // 4er ;)
        vec[i] += buffer[i];
    }
}

void accumulateVector(MeshMapping *mapping, double *vec, MPI_Comm grid) {
    int *nbors = getNbors(grid);
    accumulateVectorV(vec, grid, nbors);
    accumulateVectorE(mapping, vec, grid, nbors);
    free(nbors);
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
void accumulateVectorTrivially(MeshMapping *mapping, double *vec, MPI_Comm comm) {
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