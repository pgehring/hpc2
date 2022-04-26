#include "hpc.h"
#include <stddef.h>

index jds_spmv (const jds *A, const double *x, double *y)
{
    if (!A || !x || !y){
	return 0;
    }

    // zero initialize result vector
    for (index i=0; i<A->m; ++i){
	y[i] = 0;
    }

    index colInd;
    index rowInd;
    index nColElem = 0;

    for (index j=0; j<A->ndiag; ++j){
	index nElem = A->p[j+1] - A->p[j]; 

	for (index i=0; i<nElem; ++i){
	    colInd = A->j[nColElem+i];
	    rowInd = A->perm[i];

	    y[rowInd] += A->x[nColElem+i] * x[colInd];

	}

	nColElem += nElem; 
    }

    return (1) ;
}
