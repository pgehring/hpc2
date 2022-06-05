#ifndef _HPC_H
#define _HPC_H

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define index ptrdiff_t

/* --- primary HPC routines and data structures ------------------------- */

typedef struct sed_sparse { /* matrix in sparse matrix in compressed col. */
                            /* with extracted diagonal storage form      */
    index nzmax;            /* maximum number of entries */
    index n;                /* number of rows/columns          */
    index *i;               /* col pointers and row indices    */
    double *x;              /* numerical values, size i[n] */
} sed;

typedef struct cs_sparse /* matrix in compressed-row/col or triplet form */
{
    index nzmax; /* maximum number of entries */
    index m;     /* number of rows */
    index n;     /* number of columns */
    index *p;    /* col/row pointers (size n+1) or col indices (size nzmax) */
    index *ind;  /* row/col indices, size nzmax */
    double *x;   /* numerical values, size nzmax */
    index nz;    /* # of entries in triplet matrix,
                  * -1 for compressed-col, -2 for compressed-row */
} cs;

typedef struct mesh_data { /* mesh */
    index ncoord;          /* number of coordinates */
    index nelem;           /* number of elements */
    index nedges;          /* number of edges */
    index nbdry;           /* number of boundary elements */
    index nfixed;          /* number of fixed nodes */
    double *coord;         /* coordinates (x1,y1,x2,y2, ... ,x_ncoord,y_ncoord) */
    index *elem;           /* elements ([e1,e2,e3,m1,m2,m3,t1], ... ) */
    index *edge2no;        /* 404 */
    index *bdry;           /* bdry ([e1,e2,m1,t1], [e3,e4,m2,t2], ...) */
    index *fixed;          /* bdry ([e1,e2,m1,t1], [e3,e4,m2,t2], ...) */
} mesh;

typedef struct MeshMapping {
    mesh *localMesh;
    index *vertexL2G; /* local to global vertex mapping */
    index *elemL2G;   /* local to global elem mapping */
    index *bdryL2G;   /* local to global bdry mapping */
    index globalNcoord;
    index lMeshDimX;
    index lMeshDimY;
} MeshMapping;

// Utility functions
void *hpc_realloc(void *p, index n, size_t size, index *ok);
double hpc_cumsum(index *p, index *c, index n);

// SED functions
sed *sed_alloc(index n, index nzmax, index values);
index sed_realloc(sed *A, index nzmax);
sed *sed_free(sed *A);
sed *sed_done(sed *C, void *w, void *x, index ok);
// sed *sed_compress (const cs *A);
index sed_print(const sed *A, index brief);
index sed_gaxpy(const sed *A, const double *x, double *y);
index sed_dupl(sed *A);
index sed_gs_constr(const sed *A, const double *b, double *x, double *w,
                    index *fixed, index nFixed, index forward);

// useful CS functions for testing
index cs_entry(cs *T, index i, index j, double x);
index cs_print(const cs *A, index brief);
index cs_spmv(const cs *A, const double *x, double *y, double alpha,
              double beta);
cs *cs_alloc(index m, index n, index nzmax, index values, index typ);
cs *cs_free(cs *A);
index cs_realloc(cs *A, index nzmax);
index cs_spmv_const(const cs *A, const double *x, double *y, double alpha,
                    double beta, index *fixed, index nfixed);

// Mesh functions
mesh *mesh_alloc(index ncoord, index nelem, index nbdry);
void mesh_free(mesh *M);
mesh *mesh_load(char *fname);
index *mesh_getFixed(const index nCoord, const index *bdry, const index nBdry,
                     index *nFixed);
index mesh_print(const mesh *M, index brief);
mesh *mesh_refine(mesh *In);
index mesh_getEdge2no(const index nElem, const index *Elem, index *nEdges,
                      index **edge2no);
void mesh_getFixedNodes(MeshMapping *mapping);

// Mesh split and mapping functions
MeshMapping ***mesh_split(mesh *, int[2]);
MeshMapping *newMeshMapping(mesh *, index, index, index);
void deleteMeshMapping(MeshMapping *);
MeshMapping ***new2DMeshMapping(index, index);
void delete2DMeshMapping(MeshMapping ***, index);
MeshMapping *mesh_transfer(MeshMapping ***globalMapping, MPI_Comm grid);

// Functions for building stiffness matrix and rhs
sed *sed_sm_build(mesh *localMesh);
void mesh_build_rhs(const mesh *localMesh, double *b,
                    double (*fV)(double *, index),
                    double (*fN)(double *, index));

// Slice functions
index getSliceOffset(index, index, index);
index getSliceSize(index, index, index);
index getSliceIndex(index, index, index);

// Misc
void stima_laplace3(double p1[2], double p2[2], double p3[2], index typ,
                    double dx[6], double ax[9]);

sed *sed_nz_pattern(mesh *M);
index sed_buildS(mesh *M, sed *T);
void mesh_buildRhs(const mesh *M, double *b, double (*f)(double *, index),
                   double (*g)(double *, index));

void hpc_fmg(sed **A, double *b, double *x, index nCycle, mesh **H,
             index nLevel, index pre, index post, index gamma);
index hpc_mg(sed **A, double *b, double *x, double tol, index maxit, mesh **H,
             index nLevel, index pre, index post, index gamma);
index hpc_mg_cycle(sed **A, mesh **H, index nLevel, double **b, double **x,
                   double **r, index pre, index post, index gamma);

void hpc_rest(double *x, index *edgeno, index nEdges, double *y, index ny);
void hpc_prol(double *x, index nx, index *edgeno, index nEdges, double *y);
void hpc_prol_quad(double *x, double *y, index *elem, index nC, index nT,
                   index nE);

double kappa(double x[2], index typ);
double F_vol(double x[2], index typ);

// routines for testing
mesh *mesh_initRefinement(mesh *globalMesh, index nof_ref);

// Vector functions
double *newVector(index);
double *newVectorWithInit(index);
void accumulateVectorTrivially(MeshMapping *, double *, MPI_Comm);
void accumulateVector(MeshMapping *, double *, MPI_Comm);
void accumulateVectorV(double *, MPI_Comm, int[8]);
void accumulateVectorE(MeshMapping *, double *, MPI_Comm, int[8]);
double dot_dist(double *x_ac, double *y_dist, size_t len, MPI_Comm grid);

// Solver functions
void sed_spmv_sym(const sed *A, const double *x, double *y, double alpha,
                  double beta);
int solve_cg(MeshMapping *localMapping, sed *localSM, double *rhs, double *u_local,
              MPI_Comm grid, double tol, index maxIt);
double *accumulateResult(MeshMapping *localMapping, double *localResult,
			 MPI_Comm grid);
void insertDirichlet(double *u_local, mesh *localMesh, double (*u_D)(double *));
void blockFixedNodes(mesh *globalMesh, double *x);

#define HPC_MAX(a, b) (((a) > (b)) ? (a) : (b))
#define HPC_MIN(a, b) (((a) < (b)) ? (a) : (b))
#define HPC_CSC(A) (A && (A->nz == -1))
#define HPC_CSR(A) (A && (A->nz == -2))
#define HPC_TRIPLET(A) (A && (A->nz >= 0))

#ifdef _DEBUG
#define DEBUG_PRINT(fmt, ...)              \
    do {                                   \
        fprintf(stderr, fmt, __VA_ARGS__); \
    } while (0)
#else
#define DEBUG_PRINT(fmt, ...) (void)0
#endif

#endif
