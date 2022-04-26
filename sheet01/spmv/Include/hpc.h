#ifndef _HPC_H
#define _HPC_H
#include <assert.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>

#include <errno.h>
#include <string.h>

#include <stdbool.h>

#define index ptrdiff_t

/* --- primary HPC routines and data structures ------------------------- */

typedef struct cs_sparse /* matrix in compressed-row/col or triplet form */
{
    index nzmax ;     /* maximum number of entries */
    index m ;         /* number of rows */
    index n ;         /* number of columns */
    index *p ;        /* col/row pointers (size n+1) or col indices (size nzmax) */
    index *ind ;      /* row/col indices, size nzmax */
    double *x ;       /* numerical values, size nzmax */
    index nz ;        /* # of entries in triplet matrix, 
                       * -1 for compressed-col, -2 for compressed-row */
} cs ;

typedef struct gem_full /* general matrix form, entries stored row wise */
{
    index m ;         /* number of rows */
    index n ;         /* number of columns */
    double *x ;       /* numerical values */
} gem ;


typedef struct sky_pack  /* sym. matrix in sky storage form */
{
    index   n ;       /* number of rows/columns          */
    index  *p ;       /* col pointers (size n+1)         */
    double *d ;       /* diagonal entries (size n)       */
    double *x ;       /* off-diagonal entries, size p[n] */
} sky ;


typedef struct sed_sparse  /* matrix in sparse matrix in compressed col. */
{                          /* with extracted diagonal storage form      */
    index nzmax ;     /* maximum number of entries */
    index   n ;       /* number of rows/columns          */
    index  *i ;       /* col pointers and row indices    */
    double *x ;       /* numerical values, size i[n] */
} sed ;

typedef struct jds_sparse  /* matrix in jagged or triplet form */
{
    index m ;         /* number of rows */
    index n ;         /* number of columns */
    index ndiag ;     /* number of diagonals */
    index *p ;        /* diagonal pointers (size ndiag+1) */
    index *j ;        /* column indices, size nzmax */
    index *perm ;     /* permutation vector */
    double *x ;       /* numerical values, size nzmax */
} jds ;

/* utilities */
void *hpc_realloc (void *p, index n, size_t size, index *ok);
double hpc_cumsum (index *p, index *c, index n);
void printVector(index n, double *x);
void initVector(index n, double *x,bool random);

/* cs format */
cs *cs_alloc (index m, index n, index nzmax, index values, index typ);
index cs_realloc (cs *A, index nzmax);
cs *cs_free (cs *A);
cs *cs_done (cs *C, void *w, void *x, index ok);
index cs_spmv (const cs *A, const double *x, double *y);

cs *cs_load (FILE *f, index issym); 
index cs_entry (cs *T, index i, index j, double x);
index cs_print (const cs *A, index brief);
cs *cs_compress (const cs *T, index typ);

/* cs utils */
double randomEntry();
bool randomVoid(double lim);
void cs_RandomInit(index m, index n, double lim, cs *A);

/* gem format */
gem *gem_alloc (index n, index m);
gem *gem_free (gem *A);
index gem_gaxpy (const gem *A, const double *x, double *y);
gem *gem_compress (const cs *T);
index gem_print (const gem *A, index brief);

/* skyline format */
sky *sky_alloc (index n, index nzmax);
sky *sky_free (sky *A);
sky *sky_load (FILE *f);
index sky_print (const sky *A, index brief);
sky *sky_compress (const cs *T);
index sky_spmv (const sky *A, const double *x, double *y);

/* jagged format */
jds *jds_spalloc (index m, index n, index nz, index ndiag);
jds *jds_free (jds *A);
index jds_print (const jds *A, index brief);
jds *jds_compress (cs *A);
index jds_spmv (const jds *A, const double *x, double *y);


/* sed format */
sed *sed_alloc (index n, index nzmax, index values);
index sed_realloc (sed *A, index nzmax);
sed *sed_free (sed *A);
sed *sed_done (sed *C, void *w, void *x, index ok);
sed *sed_compress (const cs *A);
index sed_print (const sed *A, index brief);
index sed_spmv (const sed *A, const double *x, double *y);

#define HPC_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define HPC_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define HPC_CSC(A) (A && (A->nz == -1))
#define HPC_CSR(A) (A && (A->nz == -2))
#define HPC_TRIPLET(A) (A && (A->nz >= 0))
#endif
