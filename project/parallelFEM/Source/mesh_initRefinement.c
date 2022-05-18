#include "hpc.h"


/**
  * @brief perform initial refinement of the mesh to make the split possible
           (mesh dimensions must be greater than grid dimension)
  * @param globalMesh unrefined mesh
  * @param nof_ref number of refinements to be performed
*/
mesh *mesh_initRefinement(mesh *globalMesh, index nof_ref) {
    /** define initial refinement, such that the dimensions of the mesh are
        greater than the grid dimensions */
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
    for (index i = 0; i < nof_ref; ++i) {
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
