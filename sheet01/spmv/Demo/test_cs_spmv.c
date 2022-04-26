#include "hpc.h"

void print_vector(size_t n, const double *x) {
	for (size_t i = 0; i < n; i++) {
		printf("%lg\n", x[i]);
	}
}

int main() {
	cs *A = cs_load(stdin, 0);
	size_t m = A->m;
	size_t n = A->n;
	double x[n];
	double y[m];
	memset(x, 0, sizeof(double) * n);
	memset(y, 0, sizeof(double) * m);
	x[0] = 1;

    printf ("---------------------\nA:\n");
	cs_print(A, 0);

	printf ("---------------------\nx:\n");
	print_vector(A->n, x);

	assert(cs_spmv(A, x, y));
	printf ("---------------------\ny:\n");
	print_vector(A->m, y);

	return 0;
}
