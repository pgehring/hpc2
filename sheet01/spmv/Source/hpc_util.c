#include "hpc.h"

void initVector(index n, double *x,bool random){
    for (size_t i=0; i<n; ++i){
        if (random){
            x[i] = randomEntry();
        } else{
            x[i] = 0;
        }
    }
}

void printVector(index n, double *x){
    for (size_t i=0; i<n; ++i){
	    printf("%lf\n",x[i]);
    }
}