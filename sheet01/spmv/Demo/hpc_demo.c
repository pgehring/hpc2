#include "hpc.h"
int main (void)
{
    cs *T, *C, *R;
    gem *G ;
    sed *S ;
     
    T = cs_load (stdin, 0) ;             /* load triplet matrix T from stdin */
    printf ("---------------------\nT triplet:\n") ; 
    cs_print (T, 0) ;                    /* print T */
    
    C = cs_compress(T,1) ;               /* C = compressed-col form of T */     
    printf ("---------------------\nC comp. col:\n") ; 
    cs_print (C, 0) ;                    /* print C */
    
    R = cs_compress(T,2) ;               /* R = compressed-row form of T */
    printf ("---------------------\nR comp. row:\n") ; 
    cs_print (R, 0) ;                    /* print R */

    G = gem_compress(T) ;                /* G = general full matrix storage of T */
    printf ("---------------------\nG gen. full:\n") ; 
    gem_print (G, 0) ;                   /* print G */

    S = sed_compress(T) ;                /* S = sparse extr. diag. of T */
    printf ("---------------------\nS sparse diag:\n") ; 
    sed_print (S, 0) ;                  /* print S */
    
    printf ("---------------------\n\n") ; 

    cs_free (T) ;                        /* clear T */
    cs_free (C) ;                        /* clear C */
    cs_free (R) ;                        /* clear R */
    gem_free (G) ;                       /* clear G */
    sed_free (S) ;                       /* clear S */
    return (0) ;
}
