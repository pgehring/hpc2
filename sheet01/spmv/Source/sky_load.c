#include "hpc.h"
/* load a sym sky matrix from a file */
sky *sky_load (FILE *fp)
{
  double i, j, x, *Sx, *Sd;
  index ii,jj,k, offset, m, n, *Sp;
  sky *S ;
  if (!fp) return (NULL) ;                              /* check inputs */

  offset = 1024; m = 0; n = 0;
  while (fscanf (fp, "%lg %lg %lg\n", &i, &j, &x) == 3) /* get dimension */
  {                                                     /* of matrix */
    offset = HPC_MIN(offset,(index) i); m = HPC_MAX(m, (index) i);
    offset = HPC_MIN(offset,(index) j); n = HPC_MAX(n, (index) j);
  }
  if ( m != n) return (NULL) ;                          /* check dimensions */
  m = m - offset + 1 ;   n = n - offset + 1 ;
  S = sky_alloc (n, 0);
  if ( !S ) sky_free(S) ;                               /* out of memory */
  Sp = S->p; Sd = S->d;
  
  rewind(fp);                                           /* get envelope */
  for (k = 0; k < n; k++) Sp[k] = k;
  while (fscanf (fp, "%lg %lg %lg\n", &i, &j, &x) == 3) /* get envelope */
  {
    Sp[(index) i - offset] = HPC_MIN(Sp[(index) i - offset], (index) j - offset);
  }
  Sp[0] = 0;
  for (k = 1; k < n; k++) Sp[k] = Sp[k-1] + ( k - Sp[k] );

  S->x = Sx = calloc (Sp[n-1], sizeof (double)) ;
  if ( !Sx ) sky_free(S) ;                              /* out of memory */
  for (k = 0; k < n; k++) Sd[k] = 0.0;
  rewind(fp);                                           /* load data into place */
  while (fscanf (fp, "%lg %lg %lg\n", &i, &j, &x) == 3)  
  {                                                     
    ii = (index) i - offset; jj = (index) j - offset;
    if (ii == jj) Sd[ii] += x;
    else if (ii > jj) Sx[ Sp[ii]-ii+jj ] += x;
  }
  return (S) ;
}
