#include <assert.h>

#include "hpc.h"

/* allocate mesh data structure */
mesh *mesh_alloc(index ncoord, index nelem, index nbdry) {
    mesh *M = (mesh *)malloc(sizeof(mesh)); /* allocate the mesh struct */
    if (!M) {
        fprintf(stderr, "[E] Could not allocate mesh\n");
        abort();
    }

    M->ncoord = ncoord;
    M->nelem = nelem;
    M->nbdry = nbdry;
    M->nedges = 0;
    M->nfixed = 0;
    M->coord = (double *)malloc(ncoord * 2 * sizeof(double));
    M->elem = (index *)malloc(nelem * 7 * sizeof(index));
    M->bdry = (index *)malloc(nbdry * 4 * sizeof(index));
    M->edge2no = NULL;
    M->fixed = NULL;

    if (!M->coord || !M->elem || !M->bdry) {
        fprintf(stderr, "[E] Could not allocate mesh attributes\n");
        abort();
    }

    return M;
}

/* free a mesh data structure */
void mesh_free(mesh *M) {
    if (!M) return; /* do nothing if M already NULL */
    free(M->coord);
    free(M->elem);
    free(M->bdry);
    free(M->edge2no);
    free(M->fixed);
    free(M);
}

/**
 * @brief Allocates memory for a mesh mapping struct.
 *
 * @param localMesh Must not be NULL
 * @param globalNcoords Global number of vertices
 * @return MeshMapping*
 */
MeshMapping *newMeshMapping(mesh *localMesh, index globalNcoords) {
    assert(localMesh != NULL);

    MeshMapping *mapping = (MeshMapping *)malloc(sizeof(MeshMapping));
    if (!mapping) {
        fprintf(stderr, "[E] Could not allocate mapping\n");
        abort();
    }

    mapping->localMesh = localMesh;
    mapping->vertexL2G = (index *)malloc(localMesh->ncoord * sizeof(index));
    mapping->elemL2G = (index *)malloc(localMesh->nelem * sizeof(index));
    mapping->bdryL2G = (index *)malloc(localMesh->nbdry * sizeof(index));
    mapping->globalNcoord = globalNcoords;

    if (!mapping->vertexL2G || !mapping->elemL2G || !mapping->bdryL2G) {
        fprintf(stderr, "[E] Could not allocate mapping attributes\n");
        abort();
    }

    return mapping;
}

/**
 * @brief Frees the space for the mapping, but not the space for the global or
 * local meshes.
 *
 * @param mapping
 */
void deleteMeshMapping(MeshMapping *mapping) {
    if (!mapping) {
        return;
    }

    free(mapping->vertexL2G);
    free(mapping->elemL2G);
    free(mapping->bdryL2G);
    free(mapping);
}

/**
 * @brief Allocates memory for 2D-array of pointers. Does not allocate space for
 * the actual mappings.
 *
 * @param gridDimX
 * @param gridDimY
 * @return MeshMapping***
 */
MeshMapping ***new2DMeshMapping(index gridDimX, index gridDimY) {
    MeshMapping ***mapping =
        (MeshMapping ***)malloc(gridDimX * sizeof(MeshMapping **));
    if (!mapping) {
        fprintf(stderr,
                "[E] Could not allocate mapping memory (gridDimX=%zd)\n",
                gridDimX);
        abort();
    }

    for (index k = 0; k < gridDimX; k++) {
        mapping[k] = (MeshMapping **)malloc(gridDimY * sizeof(MeshMapping *));
        if (!mapping[k]) {
            fprintf(stderr,
                    "[E] Could not allocate mapping memory (gridDimY=%zd)\n",
                    gridDimY);
            abort();
        }
    }

    return mapping;
}

/**
 * @brief Frees memory of 2D-array of pointers. Does not free space for the
 * actual mappings.
 *
 * @param mapping
 * @param gridDimX
 */
void delete2DMeshMapping(MeshMapping ***mapping, index gridDimX) {
    if (!mapping) {
        return;
    }
    for (index k = 0; k < gridDimX; k++) {
        free(mapping[k]);
    }
    free(mapping);
}
