#include "hpc.h"
#include <time.h>
#include <math.h>

#ifndef MAX_DIM
#define MAX_DIM 1000 
#endif

#ifndef NUM_TEST_IT
#define NUM_TEST_IT 3
#endif

#ifndef ALPHA
#define ALPHA 1
#endif

#ifndef BETA
#define BETA 0
#endif

/**
  * @brief function to create random matrix entry
  * @return random double value
*/
double randomEntry() { return ((double)rand() - RAND_MAX / 2) / RAND_MAX *10; }

/**
  * @brief function to randomly determine wheter a matrix entry should be zero
  * @param lim probability for elements to be zero
  * @return boolean value cooresponding to the result of the decision
*/
bool randomVoid(double lim) { 
    double randVal = (double)rand() / RAND_MAX;
    bool isVoid = (randVal<=lim)?true:false;
    return isVoid;
}

/**
  * @brief function to initialize a sparse symmetric matrix in COO format 
	   with random values. Diagonal elements are all nonzero
  * @param m matrix dimension m
  * @param n matrix dimension n
  * @param lim probabilty of nondiagonal elements to be zero
  * @param A matrix in COO format
*/
void cs_RandomInitSym(index m, index n, double lim, cs *A){
    double rand;
    A->m = m;
    A->n = n;

    for (index i=0; i<m; ++i){
        for (index j=i; j<n; ++j){
            double randNumb = round(randomEntry()*100)/100; 
            if (i == j) {
                cs_entry(A, i, j, randNumb);
            } else if (!randomVoid(lim)){
                cs_entry(A, i, j, randNumb);
                cs_entry(A, j, i, randNumb);
            }
        }
    }
}

/**
  * @brief helper function to get an element of a COO matrix for given indices
	   i and j
  * @param i row index of element
  * @param j col index of element
  * @param A matrix in COO format
  * @return value of the element
*/
double getSedElement(index i, index j, cs *A){
    // Search for the right index of the data vector in the index vector
    for (index k=0; k<A->nzmax; ++k){
	// Return element when found
	if (A->ind[k] == i && A->p[k] == j){
	    return A->x[k];
	}
    }
}


/**
  * @brief function to compress a symmetric sparse matrix in COO format to 
	   SED format. Only the off-diagonal elements of the lower triangle are
	   stored.
  * @param A_Comp pointer to SED matrix struct, that the compressed content of 
	   the symetric sparse matrix should be written to
  * @param A symetric sparse matrix in COO format
*/
sed *sed_compressSym(sed *A_Comp, cs *A){
    A_Comp->n = A->m;

    // Insert diagonal elements to SED matrix
    for (index i=0; i<A->m; ++i){
	A_Comp->x[i] = getSedElement(i,i,A);
    }
    
    // Insert off diagonal elements (use symmetrie)
    index ptr = A->m+1;
    for (index j=0; j<A->n; ++j){
	A_Comp->i[j] = ptr; // write pointer to column
	
	// go through all elements
	for (index k=0; k<A->nz; ++k){
	    // if at index k of the data vector an element below the diagonal is
	    // stored, write it to the matrix
	    if (A->p[k] == j && A->ind[k] > j){
		A_Comp->i[ptr] = A->ind[k]; // write row index of element
		A_Comp->x[ptr] = A->x[k]; // write actual value
		ptr++; // update pointer
	    }
	}
    }

    A_Comp->x[A->n] = ptr-1;
    A_Comp->i[A->n] = ptr;

    return A_Comp;

}

/**
  * @brief helper function to random initialize a vector
  * @param x double pointer to memory space for the vector
  * @param n dimension of the vector
*/
void vec_randomInit(double *x, index n){
    for (index i=0; i<n; ++i){
	x[i] = randomEntry();
    }
}

/**
  * @brief helper function to print a vector
  * @param x double pointer to memory space of vector
  * @param n dimension of vector
*/
void vec_print(double *x, index n){
    for (index i=0; i<n; ++i){
	printf("%lf\n",x[i]);
    }
}

/**
  * @brief function to compute 2-Norm of the error of a test vector compared
	   to a reference
  * @param n dimension of vectors
  * @param yref reference vector
  * @param ytst test vector
  * @return double value of error norm
*/
double computeErrorNorm(index n, double *yref, double *ytst){
    double errNorm=-1;
    double resultCurrent=0;

    for (index i=0; i<n; ++i){
	resultCurrent += pow(yref[i]-ytst[i],2.0);
    }
    errNorm = sqrt(resultCurrent);
    
    return errNorm;
}




int main(){
    srand(time(NULL));

    printf("\n=== Start test_sed_spmv ===\n"); 
   
    cs *csTest;
    sed *sedTest;
    double *x, *yref, *ytst;
   
    printf("\n--- Test 1: visualizeation of results for 5x5 ---\n");
    // Allocate and initialize reference matrix in COO format
    csTest = cs_alloc(5,5,25,1,0);
    cs_RandomInitSym(5,5, 0.75, csTest);
    printf("\nReference matrix: \n");
    cs_print(csTest,0);

    // Allocate symmetric Test matrix in SED format and compress reference
    // matrix to it
    sedTest = sed_alloc(5, 25, 1);
    sed_compressSym(sedTest, csTest);
    printf("\nCompressed Test matrix: \n");
    sed_print(sedTest, 0);

    // Allocate and initialize vector x 
    x = malloc(25*sizeof(double));
    vec_randomInit(x, 5);

    // allocate result vectors
    ytst = malloc(5*sizeof(double));
    yref = malloc(5*sizeof(double));

    // Calculate spmv with reference and test implementation
    cs_spmv(csTest, x, yref, ALPHA, BETA);
    sed_spmv_sym(sedTest, x, ytst, ALPHA, BETA);
    printf("\nResult with reference implementation\n");
    vec_print(yref,5);
    printf("\nResult with test implementation\n");
    vec_print(ytst,5);

    // Compute and display Error
    double errNorm = computeErrorNorm(5, yref, ytst);
    printf("\nError: %lf\n", errNorm);
    
    // Free all used memory
    free(ytst);
    free(yref);
    free(x);
    sed_free(sedTest);
    cs_free(csTest);

    
    printf("\n--- Test 2: Error for higher dimenstions ---\n");
    printf("\n(%4s,%4s)|%7s\n","m","n","error");
    printf("-------------------\n");

    // Iteratively increase Matrix dimension (dimM,dimM) and continue ad before
    for (index dimM=100; dimM <= MAX_DIM; dimM+=100){
	csTest = cs_alloc(dimM,dimM,dimM*dimM,1,0);
	cs_RandomInitSym(dimM,dimM, 0.75, csTest); 

	sedTest = sed_alloc(dimM, dimM*dimM, 1);
	sed_compressSym(sedTest, csTest);

	x = malloc(MAX_DIM*sizeof(double));
	vec_randomInit(x, dimM);

        ytst = malloc(dimM*sizeof(double));
	yref = malloc(dimM*sizeof(double));

	cs_spmv(csTest, x, yref, ALPHA, BETA);
	sed_spmv_sym(sedTest, x, ytst, ALPHA, BETA);

	double errNorm = computeErrorNorm(dimM, yref, ytst);
	printf("(%-4td,%-4td)|%7lf\n",dimM,dimM, errNorm);

	free(ytst);
	free(yref);
	free(x);
	sed_free(sedTest);
	cs_free(csTest);
    }

    printf("\n=== End test_sed_spmv ===\n");

}
