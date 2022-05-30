#include "hpc.h"
#include <assert.h>


/**
  * @brief function to compute the sparse matrix vector product for symetric
	   matrices stored in SED format, i.e
	   y <- beta*y + alpha * A*x
  * @param A pointer to sed struct of Matrix
  * @param x double pointer to mem space of vector x
  * @param y double pointer to mem space of result vector y
  * @param alpha scalar value. No computations are performed for alpha * Ax if
	   alpha is zero
  * @param beta scalar value. beta=0 equals zero initialization of y
*/
void sed_spmv_sym(const sed *A, const double *x, double *y, double alpha,
		  double beta){
    // Check data
    assert(A && x && y);

    // Update y if beta != 1
    if (beta != 1){
	for (index j=0; j<A->n; ++j){
	    y[j] = (beta==0)?0:beta * y[j];
	}
    }

    // Perform computations only if beta != 0
    if (alpha != 0){
	/** Perform computations for diagonal elements of A first
	    -> go through A->x contiguously for all diagonal elements*/
	for (index i=0; i<A->n; ++i){
	    y[i] += A->x[i]*x[i];
	}

	/** Next, perform the computations for the off-diagonal elements of A, by
	    iteratin over the elements in a columnwise manner. Use
	    The symmetry of the data to spare compuations*/
	index rowInd, colInd;
	double aij;

	/** Iteration over the column indices
	    ->itermitted contiguous acess of first n indices of A->i */
	for (index j=0; j<A->n; ++j){
	    colInd = j; // column index of elements

	    /** Iteration over the row indices
	        ->itermitted contiguous access of memory from A->i[n+1] ongoing*/
	    for (index ptr=A->i[j]; ptr < A->i[j+1]; ++ptr){
		rowInd = A->i[ptr];
		aij = A->x[ptr];
		y[rowInd] += alpha*aij*x[colInd];
		y[colInd] += alpha*aij*x[rowInd];

	    }
	}
    }
	
    



}
