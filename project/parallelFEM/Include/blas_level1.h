#ifndef BLAS_LEVEL1_H
#define BLAS_LEVEL1_H

#include <stddef.h>
#include <stdio.h>

void blasl1_dscal(double * x, size_t len, const double alpha);

double blasl1_ddot(const double * x, const double * y, size_t len);
void blasl1_daxpy(double* x, double* y, index length, double alpha, double beta);
void blasl1_dcopy(const double* a, double* b, index length, double alpha);
void blasl1_icopy(const index* a, index* b, index length, double alpha);
void blasl1_intcopy(const int* a, int* b, index length, double alpha);
void blasl1_dprint(double * data, size_t len);
void blasl1_iprint(index* data, size_t len);
void blasl1_intprint(index* data, size_t len);

#endif //BLAS_LEVEL1_H
