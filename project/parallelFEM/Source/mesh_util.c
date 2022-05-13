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
