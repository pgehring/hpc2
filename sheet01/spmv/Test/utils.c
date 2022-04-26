#include "utils.h"

double calcErrorNorm(index m, double *y_Ref, double *y_Tst){
    double errNorm=0;
    for (index i=0; i<m; ++i){
	errNorm += pow(y_Ref[i] - y_Tst[i], 2.0);
    }
    errNorm = sqrt(errNorm);
    return errNorm;
}