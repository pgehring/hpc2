#include "hpc.h"
/* print a sparse matrix; use %g for integers to avoid differences with index */
index sky_print (const sky *A, index brief)
{
  index p, j, m, n, nzmax, nz, *Ap;
  double *Ad, *Ax ;
  
  if (!A) { printf ("(null)\n") ; return (0) ; }
  n = A->n ; Ap = A->p ; Ad = A->d ; Ax = A->x ;
    
  printf ("%g-by-%g, nnz off-diag: %g\n", (double) n, (double) n, (double) (Ap [n-1])) ;
  printf ("diagonal entries \n"); 
  for (j = 0 ; j < n ; j++)
  {
    printf ("      %g : %g\n", (double) j, Ad[j] ) ;
    if (brief && p > 10) { printf ("  ...\n") ; break ; }
  }
  printf ("off-diagonal entries (lower part) \n");  
  for (j = 1 ; j < n ; j++)
  {
    printf ("    row %g : locations %g to %g\n", (double) j, 
                (double) (Ap [j-1]), (double) (Ap [j]-1)) ;
    for (p = Ap [j-1] ; p < Ap [j] ; p++)
    {
      printf ("      %g : %g\n", (double) j+p-Ap [j] , Ax ? Ax [p] : 1) ;
      if (brief && p > 10) { printf ("  ...\n") ; return (1) ; }
    }
  }
  return (1) ;
}

