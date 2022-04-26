#include <sys/times.h>
#include <unistd.h>

#include "hpc.h"

#ifndef MIN_T
#define MIN_T   0.2
#endif

#ifndef MAX_DIM_M
#define MAX_DIM_M   80
#endif

#ifndef MAX_DIM_N
#define MAX_DIM_N   80
#endif


/* return real time in seconds since start of the process 
https://www.mathematik.uni-ulm.de/numerik/hpc/ws21/session03/page06.html
*/
double wallTime()
{
   static int ticks_per_second = 0;
   if (!ticks_per_second) {
      ticks_per_second = sysconf(_SC_CLK_TCK);
   }
   struct tms timebuf;
   /* times returns the number of real time ticks passed since start */
   return (double) times(&timebuf) / ticks_per_second;
}

int main ()
{
    // printf("---------------------\n");
    // printf("\tBENCHMARK\n");
    // printf("---------------------\n");

    cs *T_COO;
    jds *T_JDS;
    sed *T_SED;
    sky *T_SKY;
    
    double *x, *y;

    double t0;

    // allocate memory for x and y vectors and initilalize x
    x = malloc(MAX_DIM_N*sizeof(double));
    initVector(MAX_DIM_N, x, true);

    y = malloc(MAX_DIM_N*sizeof(double));

    // allocate COO Matrix
    T_COO = cs_alloc(MAX_DIM_M, MAX_DIM_N, MAX_DIM_M*MAX_DIM_N, 1, 0);

    printf("hpc_bench.c\n");
    printf("M\tN\tCOO\tJDS\tSED\tSKY\n");
    printf("===========================================\n");
    
    for (index M=10, N=10; M<=MAX_DIM_M && N<=MAX_DIM_N; M+=10, N+=10){
        // printf("\nrunning benchmark for M=%td and N=%td..\n", M, N);
        /* benchmark cs */
        int runs = 0;
        double t1 = 0;
        do {
            cs_RandomInit(M, N, 0.75, T_COO);

            double t0 = wallTime();
            cs_spmv(T_COO, x, y);
            t1 += wallTime() - t0;
            ++runs;
        } while (t1 < MIN_T);
        t1 /= runs;

        // printf("Time spent for multiplication of CS matrix\n%f s\n", t1);

        /* benchmark jds */
        runs = 0;
        double t2 = 0;
        do {
            T_JDS = jds_compress(T_COO);

            double t0 = wallTime();
            jds_spmv(T_JDS, x, y);
            t2 += wallTime() - t0;
            ++runs;
        } while (t2 < MIN_T);
        t2 /= runs;
        
        jds_free(T_JDS);
        // printf("Time spent for multiplication of JDS matrix\n%f s\n", t2);

        /* benchmark sed */
        runs = 0;
        double t3 = 0;
        do {
            T_SED = sed_compress(T_COO);

            double t0 = wallTime();
            sed_spmv(T_SED, x, y);
            t3 += wallTime() - t0;
            ++runs;
        } while (t3 < MIN_T);
        t3 /= runs;
        
        sed_free(T_SED);
        // printf("Time spent for multiplication of SED matrix\n%f s\n", t3);

        /* benchmark sky */
        // runs = 0;
        // double t4 = 0;
        // do {
        //     T_SKY = sky_compress(T_COO);

        //     double t0 = wallTime();
        //     sky_spmv(T_SKY, x, y);
        //     t4 += wallTime() - t0;
        //     ++runs;
        // } while (t4 < MIN_T);
        // t4 /= runs;

        // sky_free(T_SKY);
        // printf("Time spent for multiplication of SKY matrix\n%f s\n", t4);

        printf("%td\t%td\t%4.5f\t%4.5f\t%4.5f\n", M, N, t1, t2, t3);
        // printf("%td\t%td\t%4.5f\t%4.5f\t%4.5f\t%4.5f\n", M, N, t1, t2, t3, t4);
    }
    cs_free (T_COO);
}
