#include "hpc.h"
/* S = sym sky matrix of a sym triplet matrix T */
sky *sky_compress (const cs *T)
{
  double i, j, x, *Sx, *Sd, *Tx;
  index nz, k, m, n, *Sp, *Ti, *Tj;
  sky *S ;
  
  if (!HPC_TRIPLET (T)) return (NULL) ;          /* check inputs */
  m = T->m ; n = T->n; Tx = T->x ; Ti = T->ind ; Tj = T->p ; nz = T->nz ;
  if ( m != n ) return(NULL);
  S = sky_alloc(n, 0) ;                          /* allocate result */
  if (!S) sky_free(S) ;                          /* check inputs */
  Sp = S->p; Sd = S->d;
  for (k = 0; k < n; k++) Sp[k] = k;             /* get envelope */
  for (k = 0; k < nz; k++) Sp[ Ti[k] ] = HPC_MIN(Sp[Ti[k]], Tj[k]);   // Sp[0] = 0;
  for (k = 1; k < n; k++) Sp[k] = Sp[k-1] + ( k - Sp[k] );
  S->x = Sx = calloc (Sp[n-1], sizeof (double)) ;
  if ( !Sx ) sky_free(S) ;                       /* out of memory */
  for (k = 0; k < n; k++) Sd[k] = 0.0;
  for (k = 0; k < nz; k++)                       /* put data into place, */
  {                                              /* store strict lower part */                 
    if ( Ti[k] == Tj[k] ) Sd[Ti[k]] += Tx[k];
    else if ( Ti[k] > Tj[k] ) Sx[ Sp[Ti[k]] - Ti[k] + Tj[k] ] += Tx[k];
  }
  return (S) ;
}
