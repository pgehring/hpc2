#include "hpc.h"
gem *gem_compress (const cs *T)
{
    index j, k, m, n, nz, *Ti, *Tj;
    double *Cx, *Tx;
    gem *C ;
    if (!HPC_TRIPLET (T)) return (NULL) ;            /* check inputs */
    m = T->m ; n = T->n; Tx = T->x ; Ti = T->ind ; Tj = T->p ; nz = T->nz ;
    C = gem_alloc(m, n) ;                           /* allocate result */
    if (!C) return (NULL) ;                          /* out of memory */    
    Cx = C->x ;
    for (k = 0 ; k < nz ; k++) Cx[Ti[k]*m+Tj[k]] += Tx[k];
    return (C) ;
}
