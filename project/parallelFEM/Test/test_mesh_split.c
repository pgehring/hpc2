#include <assert.h>

#include "hpc.h"

int main() {
    printf("\n=== Start test_mesh_split ===\n");

    char fname[32] = "../Problem/problem1";
    mesh *m1 = mesh_load(fname);
    mesh_getEdge2no(m1->nelem, m1->elem, &m1->nedges, &m1->edge2no);
    m1->fixed = mesh_getFixed(m1->ncoord, m1->bdry, m1->nbdry, &m1->nfixed);
    mesh *m2 = mesh_refine(m1);
    mesh_getEdge2no(m2->nelem, m2->elem, &m2->nedges, &m2->edge2no);
    m2->fixed = mesh_getFixed(m2->ncoord, m2->bdry, m2->nbdry, &m2->nfixed);
    mesh *m3 = mesh_refine(m2);
    mesh_getEdge2no(m3->nelem, m3->elem, &m3->nedges, &m3->edge2no);
    m3->fixed = mesh_getFixed(m3->ncoord, m3->bdry, m3->nbdry, &m3->nfixed);

    index dims[2] = {3, 3};
    mesh_split(m3, dims);

    printf("Okay\n");
    printf("=== End test_mesh_split ===\n");
}
