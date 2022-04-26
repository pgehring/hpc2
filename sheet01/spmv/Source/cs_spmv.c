#include "hpc.h"
/* y = A * x + y */
index cs_spmv(const cs *A, const double *x, double *y) {
	index p, j, m, n, nz, *Ap, *Ai;
	double *Ax, tmp;

	if (!A || !x || !y) return (0); /* check inputs */

	if (HPC_CSC(A)) {
		n = A->n;
		Ap = A->p;
		Ai = A->ind;
		Ax = A->x;
		for (j = 0; j < n; j++) {
			for (p = Ap[j]; p < Ap[j + 1]; p++) {
				y[Ai[p]] += Ax[p] * x[j];
			}
		}
	} else if (HPC_CSR(A)) {
		m = A->m;
		Ap = A->p;
		Ai = A->ind;
		Ax = A->x;
		for (j = 0; j < m; j++) {
			for (p = Ap[j]; p < Ap[j + 1]; p++) {
				y[j] += Ax[p] * x[Ai[p]];
			}
		}
	} else {
		m = A->m;
		n = A->n;
		nz = A->nz;
		Ai = A->ind;
		Ap = A->p;
		Ax = A->x;
		for (j = 0; j < nz; j++) {
			y[Ai[j]] += Ax[j] * x[Ap[j]];
		}	
	}

	return (1);
}
