#include <assert.h>

#include "hpc.h"

#define GRID_M 2
#define GRID_N 2

/**
 * @brief compute the number of nodes for each local mesh
 * and write it to array
 *
 * @param buffer Array
 * @param globalMesh global Mesh
 * @param gridDimX x dimension of MPI Grid
 * @param gridDimY y dimension of MPI Grid
 */
void getTargetNodeCounts(index *buffer, mesh *globalMesh, index gridDimX,
                         index gridDimY) {
    // Compute the number of nodes for all local meshes
    for (index i = 0; i < gridDimY; ++i) {
        for (index j = 0; j < gridDimX; ++j) {
            buffer[i * gridDimX + j] =
                getSliceSize((index)sqrt((double)globalMesh->ncoord - 1),
                             gridDimY, i) +
                1;
            buffer[i * gridDimX + j] *=
                getSliceSize((index)sqrt((double)globalMesh->ncoord - 1),
                             gridDimX, j) +
                1;
        }
    }
}

/**
 * @brief Collect the numer of nodes of all local meshes and write the values to
 *	   annray
 * @param buffer array
 * @param mapping mesh mapping
 * @param gridDimX x dimension of MPI grid
 * @param gridDimY y dimension of MPI grid
 */
void collectLocalNodeCounts(index *buffer, MeshMapping ***mapping,
                            index gridDimX, index gridDimY) {
    for (int i = 0; i < gridDimY; ++i) {
        for (int j = 0; j < gridDimX; ++j) {
            buffer[i * gridDimX + j] = mapping[j][i]->localMesh->ncoord;
        }
    }
}

/**
 * @brief Compute the expected number of boundary elements for all lcoal meshes
 *	   and write the values to array
 * @param buffer array
 * @param mesh global Mesh
 * @param gridDimX x dimension of MPI grid
 * @param gridDimY y dimension of MPI grid
 */
void getTargetBndryCounts(index *buffer, mesh *mesh, index gridDimX,
                          index gridDimY) {
    index nEdgesX, nEdgesY;

    // First, compute expected number of edges for each local mesh
    for (index i = 0; i < gridDimY; ++i) {
        nEdgesY =
            getSliceSize((index)sqrt((double)mesh->ncoord - 1), gridDimY, i);

        for (index j = 0; j < gridDimX; ++j) {
            nEdgesX = getSliceSize((index)sqrt((double)mesh->ncoord - 1),
                                   gridDimX, j);

            /** determin expected number of boundary elements in dependence of
               the indices of the current local mesh.
             */
            buffer[i * gridDimX + j] = 0;

            if (i == 0 || i == gridDimY - 1) {
                buffer[i * gridDimX + j] = nEdgesX;
            }
            if (j == 0 || j == gridDimX - 1) {
                buffer[i * gridDimX + j] += nEdgesY;
            }
        }
    }
}

/**
  * @brief collect the number of boundary elements of all local meshes according
           according to the mapping and write the values to an array
  * @param buffer array
  * @param mapping mesh mapping
  * @param gridDimX x dimension of MPI grid
  * @param gridDimY y dimension of MPI grid
*/
void collectLocalBndryCounts(index *buffer, MeshMapping ***mapping,
                             index gridDimX, index gridDimY) {
    for (int i = 0; i < gridDimY; ++i) {
        for (int j = 0; j < gridDimX; ++j) {
            buffer[i * gridDimX + j] = mapping[j][i]->localMesh->nbdry;
        }
    }
}

/**
  * @brief perform initial refinement of the mesh to make the split possible
           (mesh dimensions must be greater than grid dimension)
  * @param globalMesh unrefined mesh
  * @param gridDims array of dimensions of MPI grid
*/
mesh *initRefinement(mesh *globalMesh, index gridDims[2]) {
    /** define initial refinement, such that the dimensions of the mesh are
        greater than the grid dimensions */
    index nofInitRefinements = HPC_MAX(gridDims[0], gridDims[1]);
    printf("nofInitRefinements=%zd\n", nofInitRefinements);
    mesh *mesh_current;
    mesh *mesh_previous;

    // first refinement to initialize variable mesh_current
    mesh_current = mesh_refine(globalMesh);
    mesh_getEdge2no(mesh_current->nelem, mesh_current->elem,
                    &mesh_current->nedges, &mesh_current->edge2no);
    mesh_current->fixed =
        mesh_getFixed(mesh_current->ncoord, mesh_current->bdry,
                      mesh_current->nbdry, &mesh_current->nfixed);

    // Further refinements are performed iteratively
    for (index i = 0; i < nofInitRefinements; ++i) {
        mesh_previous = mesh_current;
        // Overwrite mesh_current
        mesh_current = mesh_refine(mesh_current);
        mesh_getEdge2no(mesh_current->nelem, mesh_current->elem,
                        &mesh_current->nedges, &mesh_current->edge2no);
        mesh_current->fixed =
            mesh_getFixed(mesh_current->ncoord, mesh_current->bdry,
                          mesh_current->nbdry, &mesh_current->nfixed);
        // free the memory space of the previous mesh
        mesh_free(mesh_previous);
    }

    return mesh_current;
}

/**
  * @brief Test the partitioning of the mesh, cmoputed using the function
           mesh_split() in terms of the numer of nodes and boundary elements of
           each local mesh
  * @param globalMesh gloabl mesh
  * @param gridDims araay of dimensions of the MPI grid
  * @param numRefines number of grid refinements for wich the function shall be
            tested
*/
void testNodePartitioning(mesh *globalMesh, index *gridDims, int numRefines) {
    mesh *mesh_current;
    mesh *mesh_previous;

    // initial refinement of the mesh to make the mapping possible according to
    // the mpi grid
    mesh_current = initRefinement(globalMesh, gridDims);
    index initialRefines = HPC_MAX(gridDims[0], gridDims[1]);

    printf("\n%-4s|%-5s|%-14s|%-14s|%-8s|%-8s|%-5s\n", "ref", "(i,j)",
           "cnt nodes ref", "cnt nodes tst", "bdry ref", "bdry tst", "pass");
    printf("---------------------------------------------------------------\n");

    // Refine mesh iteratively and collect target and local node counts
    index refNodeCounts[gridDims[0] * gridDims[1]];
    index testNodeCounts[gridDims[0] * gridDims[1]];
    index refBndryCounts[gridDims[0] * gridDims[1]];
    index tstBndryCounts[gridDims[0] * gridDims[1]];
    MeshMapping ***testMapping;

    for (index i = initialRefines; i < numRefines; ++i) {
        // create mapping
        testMapping = mesh_split(mesh_current, gridDims);

        // Get raference node counts and node counts of local meshes
        getTargetNodeCounts(refNodeCounts, mesh_current, gridDims[0],
                            gridDims[1]);
        collectLocalNodeCounts(testNodeCounts, testMapping, gridDims[0],
                               gridDims[1]);

        // Get refenerence and local bndry edge counts
        getTargetBndryCounts(refBndryCounts, mesh_current, gridDims[0],
                             gridDims[1]);
        collectLocalBndryCounts(tstBndryCounts, testMapping, gridDims[0],
                                gridDims[1]);

        // Compare and print result
        index refCount, tstCount, refBndryCnts, tstBndryCnts;
        for (index k = 0; k < gridDims[1]; ++k) {
            for (index l = 0; l < gridDims[0]; ++l) {
                refCount = refNodeCounts[k * gridDims[0] + l];
                tstCount = testNodeCounts[k * gridDims[0] + l];
                refBndryCnts = refBndryCounts[k * gridDims[0] + l];
                tstBndryCnts = tstBndryCounts[k * gridDims[0] + l];
                printf("%-4td|(%td,%td)|%-14td|%-14td|%-8td|%-8td|%-5s\n", i, k,
                       l, refCount, tstCount, refBndryCnts, tstBndryCnts,
                       (refCount == tstCount && refBndryCnts == tstBndryCnts)
                           ? "true"
                           : "false");
            }
        }
        printf(
            "---------------------------------------------------------------"
            "\n");
        // next refinement
        mesh_previous = mesh_current;
        mesh_current = mesh_refine(mesh_previous);
        mesh_getEdge2no(mesh_current->nelem, mesh_current->elem,
                        &mesh_current->nedges, &mesh_current->edge2no);
        mesh_current->fixed =
            mesh_getFixed(mesh_current->ncoord, mesh_current->bdry,
                          mesh_current->nbdry, &mesh_current->nfixed);

        // delete old mapping
        delete2DMeshMapping(testMapping, gridDims[0]);

        // free old mesh
        mesh_free(mesh_previous);
    }

    // free final mesh
    mesh_free(mesh_current);
}

int main() {
    setbuf(stdout, NULL);  // Disable stdout buffering
    printf("\n=== Start test_mesh_split ===\n");

    char fname[32] = "../Problem/problem1";
    mesh *m1 = mesh_load(fname);
    mesh_getEdge2no(m1->nelem, m1->elem, &m1->nedges, &m1->edge2no);
    m1->fixed = mesh_getFixed(m1->ncoord, m1->bdry, m1->nbdry, &m1->nfixed);

    printf("\n--- Test 1: Partitioning of Nodes ---\n");
    index dims[2] = {2, 4};

    testNodePartitioning(m1, dims, 8);

    printf("Okay\n");
    printf("=== End test_mesh_split ===\n");

    mesh_free(m1);
}
