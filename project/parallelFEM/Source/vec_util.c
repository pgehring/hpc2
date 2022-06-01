#include "hpc.h"

/**
 * @brief Initializes the vector.
 *
 * @param n Length of vec
 * @param vec Pointer to vector
 * @param init Initial value
 */
void vectorInit(index n, double *vec, double init) {
    for (index i = 0; i < n; i++) {
        vec[i] = init;
    }
}

/**
 * @brief Allocates memory for a double array of size n.
 * Aborts if there is not sufficient memory.
 *
 * @param n Size of the allocated array
 * @return double* Pointer to array
 */
double *newVector(index n) {
    double *vec = (double *)calloc(n, sizeof(double));
    if (!vec) {
        printf("[E] Could not allocate memory for vector\n");
        abort();
    }
    return vec;
}

/**
 * @brief Allocates memory for a double array of size n
 * and initializes it to 0.
 * Aborts if there is not sufficient memory.
 *
 * @param n Size of the allocated array
 * @return double* Pointer to array
 */
double *newVectorWithInit(index n) {
    double *vec = newVector(n);
    vectorInit(n, vec, 0);
    return vec;
}
