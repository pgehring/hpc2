#include "hpc.h"
/* allocate a sparse matrix in jagged diagonal storage format */
jds *jds_spalloc (index m, index n, index nz, index ndiag)
{
    jds *A = calloc (1, sizeof (jds)) ;      /* allocate the jds struct */
    if (!A) return (NULL) ;                  /* out of memory */
    A->m = m ;                               /* define dimensions and nzmax */
    A->n = n ;
    A->ndiag = ndiag ;
    A->p = malloc ((ndiag+1)*sizeof (index)) ;
    if (A->p) A->p[ndiag] = nz = HPC_MAX (nz, 1);
    A->j = malloc (nz * sizeof (index)) ;
    A->x = malloc (nz * sizeof (double));
    A->perm = malloc (m * sizeof (index)) ;
    return ((!A->j || !A->x || !A->perm ) ? jds_free (A) : A) ;
}

jds *jds_free (jds *A)
{
    if (!A) return (NULL) ;     /* do nothing if A already NULL */
    free (A->p) ;
    free (A->j) ;
    free (A->perm) ;
    free (A->x) ;
    free (A) ;
    return (NULL) ;   /* free the jds struct and return NULL */
}



