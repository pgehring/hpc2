#include <assert.h>
#include <mpi.h>

#include "hpc.h"

void printTstArray(index n, double* vec) {
    fprintf(stderr, "tst:");
    for (int i = 0; i < n; i++) {
        fprintf(stderr, "%12.0lf ", vec[i]);
    }
    fprintf(stderr, "\n");
}

void printRefArray(index n, double* vec, index* vertexL2G) {
    fprintf(stderr, "ref:");
    for (int i = 0; i < n; i++) {
        fprintf(stderr, "%12.0lf ", vec[vertexL2G[i]]);
    }
    fprintf(stderr, "\n");
}

void testAccumulateVector(MPI_Comm grid) {
    int dims[2];
    int periods[2];
    int coords[2];
    int rank;
    int nofProcesses;
    MPI_Cart_get(grid, 2, dims, periods, coords);
    MPI_Comm_rank(grid, &rank);
    MPI_Comm_size(grid, &nofProcesses);

    // Load mesh
    char fname[32] = "../Problem/problem1";
    mesh* m1 = mesh_load(fname);
    mesh_getEdge2no(m1->nelem, m1->elem, &m1->nedges, &m1->edge2no);
    m1->fixed = mesh_getFixed(m1->ncoord, m1->bdry, m1->nbdry, &m1->nfixed);

    // Refine mesh
    index nofRefinements = 7;
    mesh* m2 = mesh_initRefinement(m1, nofRefinements);
    MeshMapping*** mappings = mesh_split(m2, dims);

    // Local mesh and mappings
    MeshMapping* localMapping = mappings[coords[0]][coords[1]];
    mesh* localMesh = localMapping->localMesh;
    index gncoord = localMapping->globalNcoord;
    index lncoord = localMesh->ncoord;

    // Construct vector to accumulate
    double* vec = newVector(lncoord);
    for (index i = 0; i < lncoord; i++) {
        vec[i] = rank * 1e+3 + i;
    }

    accumulateVector(localMapping, vec, grid);
    // accumulateVectorTrivially(localMapping, vec, grid);

    // Construct reference vector
    double* ref = newVectorWithInit(gncoord);
    for (index r = 0; r < nofProcesses; r++) {
        MPI_Cart_coords(grid, r, 2, coords);
        MeshMapping* mapping = mappings[coords[0]][coords[1]];
        for (index i = 0; i < mapping->localMesh->ncoord; i++) {
            ref[mapping->vertexL2G[i]] += r * 1e+3 + i;
        }
    }

    // Assert that it is the same
    for (index i = 0; i < lncoord; i++) {
        assert(vec[i] == ref[localMapping->vertexL2G[i]]);
    }

    free(ref);
    free(vec);
    mesh_free(m2);
    mesh_free(m1);
}

// void testAccumulateVectorV(MPI_Comm grid) {
//     int dims[2];
//     int periods[2];
//     int coords[2];
//     MPI_Cart_get(grid, 2, dims, periods, coords);

//     double vec[4] = {1, 10, 100, 1000};
//     accumulateVectorV(vec, grid);

//     fprintf(stderr, "(%d,%d)\n", coords[0], coords[1]);
//     fprintf(stderr, "%6.0lf - %6.0lf\n", vec[2], vec[3]);
//     fprintf(stderr, "%6.0lf - %6.0lf\n", vec[0], vec[1]);
// }

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    // Set up 2-dim grid
    int nofProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &nofProcesses);
    int dims[2] = {0, 0};
    MPI_Dims_create(nofProcesses, 2, dims);
    int periods[2] = {0, 0};
    MPI_Comm grid;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid);
    int rank;
    MPI_Comm_rank(grid, &rank);
    int coords[2];
    MPI_Cart_coords(grid, rank, 2, coords);
    if (rank == 0) {
        printf("\n=== Start test_vec_accumulate ===\n");
    }

    testAccumulateVector(grid);

    MPI_Finalize();
    if (rank == 0) {
        printf("--- Total Result: PASS ---\n");
        printf("=== End test_vec_accumulate ===\n");
    }
}
