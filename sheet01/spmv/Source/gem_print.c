#include "hpc.h"
/* print a full matrix; use %g for integers to avoid differences with index */
index gem_print (const gem *A, index brief)
{
  index i, j, m, n, nzmax, nz, *Ap, *Ai ;
  double *Ad, *Ax ;
  
  if (!A) { printf ("(null)\n") ; return (0) ; }
  n = A->n ; m = A->m ; Ax = A->x ;
    
  printf ("%g-by-%g\n", (double) n, (double) m) ;
  printf ("matrix entries row wise \n");  
  for (i = 0 ; i < n ; i++)
  {
    for (j = 0 ; j < m ; j++)
    {
      printf (" %5.3g", Ax[i*m+j]) ;
      if (brief && j > 10) { printf ("  ...\n") ; return (1) ; }
    }
    if (brief && i > 10) { printf ("  ...\n") ; return (1) ; }
    printf ("\n") ;
  }
  return (1) ;
}

