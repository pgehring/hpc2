#include "hpc.h"
/* allocate a general (n x m) matrix, entries are stored row wise */
gem *gem_alloc (index n, index m)
{
    gem *A = calloc (1, sizeof (gem)) ;    /* allocate the sky struct */
    if (!A) return (NULL) ;                /* out of memory */
    A->n = n ;                             /* define dimensions */
    A->m = m ;                             
    A->x = calloc (n * m, sizeof (double)) ; /* allocate matrix */
    return ( !(A->x) ? gem_free (A) : A) ;
}

/* free a gem matrix */
gem *gem_free (gem *A)
{
    if (!A) return (NULL) ;      /* do nothing if A already NULL */
    free (A->x) ;                /* free the gem struct and return NULL */
    free (A);
    return (NULL) ; 
}

