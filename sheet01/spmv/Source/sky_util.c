#include "hpc.h"
/* allocate a sparse matrix (triplet form or compressed-col/row form) */
sky *sky_alloc (index n, index nzmax)
{
    sky *A = calloc (1, sizeof (sky)) ;    /* allocate the sky struct */
    if (!A) return (NULL) ;                /* out of memory */
    A->n = n ;                             /* define dimensions and nzmax */
    A->p = malloc ( n * sizeof (index)) ;
    A->p[n-1] = nzmax ;
    A->d = malloc (n * sizeof (double)) ;
    A->x = nzmax ? malloc (nzmax * sizeof (double)) : NULL ;
    return ((!A->p ||!A->d || (nzmax && !A->x)) ? sky_free (A) : A) ;
}

/* free a sparse matrix */
sky *sky_free (sky *A)
{
    if (!A) return (NULL) ;      /* do nothing if A already NULL */
    free (A->p) ;                /* free the crs struct and return NULL */
    free (A->x) ;
    free (A);
    return (NULL) ; 
}

/* free workspace and return a skyline matrix result */
sky *sky_done (sky *C, void *w, void *x, index ok)
{
    free (w) ;                         /* free workspace */
    free (x) ;
    return (ok ? C : sky_free (C)) ;   /* return result if OK, else free it */
}


