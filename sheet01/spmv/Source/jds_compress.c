#include "hpc.h"
/* computes a permutation (by the distribution counting sort) which, when 
 * applied to the input vector v, sorts the integers in v in descending order 
 * Reference: Donald Knuth, The Art of Computer Programming, 
 * Volume 3: Sorting and Searching, Addison-Wesley, 1973, pages 78-79. */
index dcsort (index *v, index *p, index n, index *w, index low, index high)
{
  index i, vi;
  /* dim:  p[n], w[high-low+1] */
  for (i = 0 ; i < high-low+1 ; i++) w[i] = 0;
  for (i = 0 ; i < n ; i++) w[v[i]-low]++;
  w[high-low]--;  
  for (i = high-low ; i >0 ; i--) w[i-1]+=w[i];
  for (i = n-1 ; i >=0 ; i--){
    vi = v[i]-low;
    p[w[vi]] = i;
    w[vi]--;
  }
  return (1);
}

/* converts compressed sparse row to jagged diagonal storage */
jds *jds_compress(cs *T)
{
    index i, j, jj, k, m, n, nz, *Ap, *Aj ;
    index low, high, len, *Cp, *Cperm, *Cj, ndiag=0, *Ti;
    jds *C;
    cs *A;
    double *Ax, *Cx ;
    
    A = HPC_TRIPLET(T) ? cs_compress(T,2) : T;
    if (!HPC_CSR (A)) {cs_free(A); return (NULL) ;}  
    m = A->m ; n = A->n ; Ap = A->p ; Aj = A->ind ; Ax = A->x ; nz = Ap[m];
    /* Define initial permutation P and get lengths of each row. */
    Ti = malloc (m * sizeof (index));
    low = n;
    for (i = 0 ; i < m ; i++)
    {
      len = Ap[i+1] - Ap[i];
      ndiag = HPC_MAX( ndiag, len );
      low   = HPC_MIN( low, len );
      Ti[i] = len;
    }
    C = jds_spalloc (m, n, nz, ndiag);
    Cp = C->p ; Cperm = C->perm ; Cj = C->j ; Cx = C->x ;
    /* sort to get permutation */
    dcsort(Ti,Cperm,m,Cp,low,ndiag);
    /* Compute lengths of the jagged diagonals. */
    for (k = 0; k < ndiag; k++) Cp[k] = 0;  
    for (k = 0; k < m; k++)   
    {
      len = Ti[Cperm[k]];                
      for (i = 0; i < len; i++) Cp[i]++;
    } 
    /* Get the output matrix itself. */
    j = 0;
    for (jj = 0; jj < ndiag; jj++){
      len = Cp[jj];
      Cp[jj] = j;
      for (k = 0; k < len; k++){
        i = Ap[Cperm[k]] + jj;
        Cx[j] = Ax[i]; Cj[j] = Aj[i];
        j++;
      }
    }
    Cp[ndiag] = j;
  
    free(Ti);
    cs_free(A);
    return (C);
}
