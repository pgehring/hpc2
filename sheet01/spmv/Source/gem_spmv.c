#include "hpc.h"
/* general (n x m) matrix, entries stored row wise, y = y + A*x  */
index gem_gaxpy (const gem *A, const double *x, double *y)
{
    index j, k, m, n ;
    double *Ax;
    if (!x || !y) return (0) ;          /* check inputs */
    m = A->m ; n = A->n ; Ax = A->x ;
    for (j = 0 ; j < n ; j++,Ax += m)
    {
      for (k = 0 ; k < m ; k++) y[j] += Ax[k] * x[k];
    }
    return (1) ;
}
