#include "hpc.h"
/* print a matrix in jagged diagonal storage; 
 * use %g for integers to avoid differences with index */
index jds_print (const jds *A, index brief)
{
    index p, i, j, m, n, ndiag, *Ap, *Aj, *Aperm ;
    double *Ax ;
    if (!A) { printf ("(null)\n") ; return (0) ; }
    m = A->m ; n = A->n ; Ap = A->p ; Aj = A->j ; Ax = A->x ; Aperm = A->perm ;
    ndiag = A->ndiag;
  
    printf ("%g-by-%g, nnz: %g\n", (double) m, (double) n, (double) (Ap[ndiag]) ) ;
    for (j = 0 ; j < ndiag ; j++)
    {
      printf ("   diag %g : locations %g to %g\n", (double) j, 
                 (double) (Ap[j]), (double) (Ap[j+1]-1)) ;
      for (p = Ap[j], i=0 ; p < Ap[j+1] ; p++, i++)
      {
        printf ("     (%g, %g) : %g\n", (double) (Aperm[i]), 
                                   (double) (Aj [p]), Ax [p]) ;
          if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
      }
    }
    return (1) ;
}
