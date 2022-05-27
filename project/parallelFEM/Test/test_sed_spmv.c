#include "hpc.h"
#include <time.h>

#ifndef M_DIM
#define M_DIM 5 
#endif



double randomEntry() { return ((double)rand() - RAND_MAX / 2) / RAND_MAX *10; }

bool randomVoid(double lim) { 
    double randVal = (double)rand() / RAND_MAX;
    bool isVoid = (randVal<=lim)?true:false;
    return isVoid;
}

void cs_RandomInitSym(index m, index n, double lim, cs *A){
    double rand;
    
    for (index i=0; i<m; ++i){
        for (index j=i; j<n; ++j){
            double randNumb = randomEntry(); 
            if (i == j) {
                cs_entry(A, i, j, randNumb);
            } else if (!randomVoid(lim)){
                cs_entry(A, i, j, randNumb);
                cs_entry(A, j, i, randNumb);
            }
        }
    }
}

index getUniqueRandomInt(index range, bool *flags){
    index ind, rowInd;

    // Search for index in vector of row indices that is not already taken 
    while(true){
	ind = rand()%range;
	if (!flags[ind]){
	    flags[ind] = 1;
	    break;
	}
    }

    return ind;

}

double getSedElement(index i, index j, cs *A){
    // Search for the right index of the data vector in the index vector
    for (index k=0; k<A->nzmax; ++k){
	// Return element when found
	if (A->ind[k] == i && A->p[k] == j){
	    return A->x[k];
	}
    }
}


sed *sed_compressSym(cs *A){
    sed *A_Comp = sed_alloc(A->m, A->nzmax, 1);
    
    // Insert diagonal elements to SED matrix
    for (index i=0; i<A->m; ++i){
	A_Comp->x[i] = getSedElement(i,i,A);
    }
    
    // Insert off diagonal elements (use symmetrie)
    index ptr = A->m+1;
    for (index j=0; j<A->n; ++j){
	A_Comp->i[j] = ptr; // write pointer to column
	
	// go through all elements
	for (index k=0; k<A->nzmax; ++k){
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


int main(){
    srand(time(NULL));
    // Allocate and initialize matrix in triplet (COO) format for testing
    cs *csTest = cs_alloc(M_DIM,M_DIM,M_DIM*M_DIM,1,0);
    cs_RandomInitSym(M_DIM,M_DIM, 0.75, csTest); 
    
    printf("\nRef Matrix in COO-Format\n");
    cs_print(csTest,0);

    // Compress to symmetric sed format
    sed *sedTest = sed_compressSym(csTest);
    
    printf("\nSymmetric commpressed Test Matrix in SED Format\n");
    sed_print(sedTest,0);

    sed_free(sedTest);
    cs_free(csTest);
}
