#include "hpc.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef M
#define M 10
#endif
#ifndef N
#define N 10
#endif
#ifndef LIM
#define LIM 0.75
#endif

int main(){
    cs *T_COO;
    jds *T_JDS;
    gem *T_GEM;
    sky *T_SKY;

    double *x;
    double *y_Ref;
    double *y_Tst;

    printf ("---------------------\nTest jds_spmv:\n") ; 

    // initialize random value engine
    srand(time(NULL));

    // Init Test Matrix in COO format with random values
    T_COO = cs_alloc(M, N, M*N, 1, 0);
    cs_RandomInit(M, N, LIM, T_COO);

    // Create Matrix in general form for reference implementation
    T_GEM = gem_compress(T_COO);
    printf("\nTest Matrix in general Form\n");
    gem_print(T_GEM,0);

    // Create Matrix in JDS format
    T_JDS = jds_compress(T_COO);
    printf("\nTest Matrix JDS-Format:\n");
    jds_print(T_JDS, 0);

    // Create Matrix in SKY format
    T_SKY = sky_compress(T_COO);
    printf("\nTest Matrix SKY-Format:\n");
    sky_print(T_SKY, 0);

    // allocate and random initialize Vector x
    x = malloc(N*sizeof(double));
    initVector(N, x, true);
    printf("\nTest Vector:\n");
    printVector(N, x);

    // allocate and zero initialize reference result vector
    y_Ref = malloc(M*sizeof(double));
    initVector(M, y_Ref, false);
    
    // allocate and zero initialize test result vector
    y_Tst = malloc(M*sizeof(double));
    initVector(M, y_Tst, false);

    // Copmute with using reference implementation
    gem_gaxpy(T_GEM, x, y_Ref);
    printf("\nReference Result:\n");
    printVector(M, y_Ref);
    
    // Compute using jds implementation
    jds_spmv(T_JDS, x, y_Tst);    
    printf("\nTest Result:\n");
    printVector(M, y_Tst);

    //Calculate RMS error
    double err = calcErrorNorm(M, y_Ref, y_Tst);
    printf("\nRMS Error: %lf\n", err);

    // free memory
    free(y_Tst);
    free(y_Ref);
    free(x);
    cs_free(T_COO);
    gem_free(T_GEM);
    jds_free(T_JDS);

    printf("done.\n");
}
