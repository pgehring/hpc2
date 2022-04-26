#include "hpc.h"

double drand() { return ((double)rand() - RAND_MAX / 2) / RAND_MAX; }

void print_vector(size_t n, const double *x) {
	for (size_t i = 0; i < n; i++) {
		printf("%lg\n", x[i]);
	}
}

int main() {
	size_t m = 10;
	size_t n = 10;
	gem *G = gem_alloc(m, n);
	cs *C = cs_alloc(m, n, m * n, m * n, 0);
	double *x = calloc(n, sizeof(double));
	double *yref = calloc(m, sizeof(double));
	double *ytst = calloc(m, sizeof(double));

	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < n; j++) {
			G->x[i * n + j] = drand();
			cs_entry(C, i, j, G->x[i * n + j]);
		}
	}
	for (size_t j = 0; j < n; j++) {
		x[j] = drand();
	}
	memset(yref, 0, sizeof(double) * m);
	memset(ytst, 0, sizeof(double) * m);
    printf ("---------------------\nTest cs_spmv:\n") ; 
	if (n <= 10) gem_print(G, 0);

	size_t nz = C->nz;
	index *rows = C->ind;
	index *cols = C->p;
	for (size_t i = 0; i < nz; i++) {
		assert(G->x[rows[i] * n + cols[i]] == C->x[i]);
	}

	gem_gaxpy(G, x, yref);
	cs_spmv(C, x, ytst);
	for (size_t i = 0; i < m; i++) {
		assert(yref[i] == ytst[i]);
	}
	printf("x:\n");
	print_vector(n, x);
	printf("y:\n");
	print_vector(m, yref);
	printf("Done.\n");

	return 0;
}
