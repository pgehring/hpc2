#include "hpc.h"
/* solve Ax=b for x, using forward/backward substitution.
   A has L*U representation and stored in gem format */
index gem_gausssol(gem *A, double *x) 
{
  index   i, j, n ;
  double  *Ax ;
  if (!A) { printf ("(null)\n") ; return (0) ; }
  n = A->n ; Ax = A->x ; 
  for( i = 0; i < n; i++, Ax+=n )           /* forward loop */
  {
    for (j = 0 ; j < i; j++) x[i] -= Ax[j] * x[j];
  }
  for( i = n-1, Ax-=n; i >= 0; i--, Ax-=n )        /* backward loop */
  {
    for (j= i+1; j < n; j++) x[i] -= Ax[j] * x[j];
    x[i]/=Ax[i];
  }
  return(1);
}
