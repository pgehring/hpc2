#include <mpi.h>

#include "hpc.h"

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