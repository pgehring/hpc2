#include "hpc.h"
/* perform Gauss decomposition without pivoting on matrix stored 
 * in general matrix format (row wise), the original matrix destroyed */

index gem_gauss(gem *A) 
{
  index   i, j, k, n ;
  double  recip, *Ax ;
  if ( (!A) || (A->m != A->n)) { printf ("(null)\n") ; return (0) ; }
  n = A->n; Ax = A->x; 
  /* consider k-th row */
  for (k = 0; k < n-1 ; k++){
    /* modify entries in A, i.e. create LU decomp*/
    recip = 1.0/Ax[k*n+k];
    for (i = k+1; i < n; i++){ 
      /* compute  l_ik, store in A */
      Ax[i*n+k] *= recip;  
      for (j = k+1; j < n; j++){
        /* mofify (n-k-1)*(n-k-1) submatrix */
        Ax[i*n+j] -= Ax[i*n+k] * Ax[k*n+j];
      }
    }
  }
  return(1);
}
