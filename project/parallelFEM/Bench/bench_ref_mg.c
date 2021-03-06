#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include "hpc.h"

#define TIME_SAVE(j,tv) (gettimeofday(&tv[j], (struct timezone *)0))
#define TIME_ELAPSED(j, k, tv) \
    (1.E+6 * (tv[k].tv_sec - tv[j].tv_sec) + (tv[k].tv_usec - tv[j].tv_usec))

double kappa(double x[2], index typ)
{
    return (1.0);
}

double F_vol(double x[2], index typ) { return (0.0); }

double g_Neu(double x[2], index typ) { return (x[0] * x[1]); }

double u_D(double x[2])
{
    //  return ( 0.0 );
    return (x[0] * x[1]);
}

void saveBenchResults(struct timeval *tv, int n, index nofDoF, int nIter, FILE *benchFile){
    double durationVals[n-1];
    for (int i=1; i<n; ++i){
	durationVals[i-1] = 1.E+3*(tv[i].tv_sec - tv[0].tv_sec)+1.E-3*(tv[i].tv_usec -
		       tv[0].tv_usec);
    }
    
    fprintf(benchFile,"%10td ",nofDoF); 
    fprintf(benchFile, "%10d ",nIter);

    for (int i=0; i<n-2; ++i){
	printf("%lf\n",durationVals[i]);
	fprintf(benchFile, "%10lf ",durationVals[i]);
    }

    fprintf(benchFile, "%10lf\n", durationVals[n-2]);
}



int benchPoissonRefMG(char *fname, int numRefines, double(*fV)(double *,index),
			      double(*fN)(double *, index), struct timeval *tv){
    
    index n, k, ncoord, nelem, nbdry, nfixed, nedges, total, *bdry,
        cnt = 0, N = 0, MAX_ITER = 100;
    double *b, *x, *w, *Coord, x1[2], x2[2], m[2];
    mesh **H, *T;
    sed **A;

    N = numRefines;

    // Save initial timepoint
    TIME_SAVE(0,tv);
    
    /* Allocate memory for hierachy */
    H = (mesh **)malloc((N + 1) * sizeof(mesh *));
    A = (sed **)malloc((N + 1) * sizeof(sed *));

    /* Load problem */
    H[0] = mesh_load(fname); /* load geometry */
    mesh_getEdge2no(H[0]->nelem, H[0]->elem, &H[0]->nedges, &H[0]->edge2no);
    H[0]->fixed =
        mesh_getFixed(H[0]->ncoord, H[0]->bdry, H[0]->nbdry, &H[0]->nfixed);

    /* Build stiffness matrix, refine mesh and create hierachy  */
    k = 0;
    while (1)
    {
        A[k] = sed_nz_pattern(H[k]); /* get pattern of matrix */
        if (!A[k]){
            return -1;
	} 
        if (!sed_buildS(H[k], A[k])){
            return -1; /* assemble coefficient matrix */
	}
	if (k >= N){
            break;
	}	    
        H[k + 1] = mesh_refine(H[k]);
        mesh_getEdge2no(H[k + 1]->nelem, H[k + 1]->elem, &H[k + 1]->nedges,
                        &H[k + 1]->edge2no);
        H[k + 1]->fixed = mesh_getFixed(H[k + 1]->ncoord, H[k + 1]->bdry,
                                        H[k + 1]->nbdry, &H[k + 1]->nfixed);
        k++;
    }
  

    n = A[N]->n;
    x = (double *)calloc(n, sizeof(double)); /* get workspace for sol*/
    w = (double *)calloc(n, sizeof(double)); /* get temporary workspace */
    b = (double *)calloc(n, sizeof(double)); /* get workspace for rhs*/
    mesh_buildRhs(H[N], b, F_vol,
                  g_Neu); /* build rhs (volume and Neumann data */
   

    /* incorporate Dirichlet data */
    ncoord = H[N]->ncoord;
    nelem = H[N]->nelem;
    nbdry = H[N]->nbdry;
    bdry = H[N]->bdry;
    H[N]->fixed = mesh_getFixed(ncoord, bdry, nbdry, &H[N]->nfixed);
    nfixed = H[N]->nfixed;
    nedges = H[N]->nedges;
    Coord = H[N]->coord;

    for (k = 0; k < nbdry; k++)
    {
        if (!bdry[4 * k + 3])
        {
            x1[0] = Coord[2 * bdry[4 * k]];
            x1[1] = Coord[2 * bdry[4 * k] + 1];
            x2[0] = Coord[2 * bdry[4 * k + 1]];
            x2[1] = Coord[2 * bdry[4 * k + 1] + 1];
            m[0] = (x1[0] + x2[0]) / 2.0;
            m[1] = (x1[1] + x2[1]) / 2.0;
            x[bdry[4 * k]] = u_D(x1);
            x[bdry[4 * k + 1]] = u_D(x2);
            x[ncoord + bdry[4 * k + 2]] = u_D(m);
        }
    }

    // Save timepoint after setting up problem data and before solving
    TIME_SAVE(1,tv);
    cnt = hpc_mg(A, b, x, 1e-10, 50, H, N, 2, 2, 1);
    // Save time after solving
    TIME_SAVE(2,tv);

    for (k = 0; k < HPC_MIN(10, A[N]->n); k++)
    {
        printf(" x[%g] = %g\n", (double)k, x[k]);
    }


    for (k = 0; k <= N; k++)
    {
        mesh_free(H[k]);
        sed_free(A[k]);
    }
    free(H);
    free(A);

    return cnt;

}

index getNumDegreesOfFreedom(int nrefines){
    index nNodes=2;
    index nDof=0;
    for (int i=2; i<=nrefines; ++i){
	nNodes*=2;
    }
    nNodes+=1;
    nDof=nNodes*nNodes;
    nDof+=2*(nNodes-1)*(nNodes)+(nNodes-1)*(nNodes-1);

    return nDof;

}

int main(int argc, char *argv[]){
    char fname_p1[32] = "../Problem/problem1";
    struct timeval tv[5];
    FILE *bench_results = fopen("bench_results/bench_results_ref_mg.csv","w+");
    int refineMax = 2;
    index nDOF;
    int nIter;

    if (argc > 1){
	refineMax = atoi(argv[1]);
    }

    fprintf(bench_results, "%10s %10s %10s %10s\n",
	    "DOFs","Iterations","T_Setup","T_Res");
    for (int nrefine=1; nrefine<=refineMax; ++nrefine){
	printf("Refinement %d\n",nrefine);

	nIter = benchPoissonRefMG(fname_p1, nrefine, F_vol, g_Neu, tv);
	
	nDOF = getNumDegreesOfFreedom(nrefine);
	printf("nDOF=%td\n",nDOF);
	printf("nIter=%d\n",nIter);

	saveBenchResults(tv, 3, nDOF, nIter, bench_results);
    }
}
