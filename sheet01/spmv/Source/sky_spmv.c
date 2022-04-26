#include "hpc.h"

/* y = A * x + y */
index sky_spmv (const sky *A, const double *x, double *y)
{
   index p, j, m, n, *Ap;
   double *Ad, *Ax;

   if (!A || !x || !y){
      return(0);
   }

   n = A->n;
   Ap = A->p;
   Ad = A->d;
   Ax = A->x;

   for (j=0; j<n; j++) {
      y[j] += Ad[j]*x[j];
   }
   for (j=1; j<n; j++){
      for (p=Ap[j-1]; p<Ap[j]; p++){
         y[j] += Ax[p]*x[j+p-Ap[j]];
	    // Da symmetrisch:
	     y[j+p-Ap[j]] += Ax[p]*x[j];
      }
   }   

   return (1) ;
}
