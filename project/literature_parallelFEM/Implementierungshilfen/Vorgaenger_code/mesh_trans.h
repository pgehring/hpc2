#ifndef _MESH_TRANS_H
#define _MESH_TRANS_H

#include "hpc.h"
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>

/**
 * @brief Kommunikationsklasse, mit der das Mesh upgedated wird.
 * @file mesh_trans.h
 * 
 */

//---------------------------------KLASSE---------------------------------------
/**
 * @brief Enthält Daten, die für die Kommunikation zwischen den einzelnen Domains
 *        notwendig ist, sowie das lokale mesh.
 * 
 */
typedef struct mesh_transfer_data {
    index ncoord_loc;                /** @brief Anzahl Koordinaten lokal */
    index ncoord_glo;                /** @brief Anzahl Koordinaten global */
    index nelem_loc;                 /** @brief Anzahl Elemente lokal */
    index nbdry_loc;                 /** @brief Anzahl Gebietsrandpunkte im Gebiet */

    double* domcoord;                /** @brief lokale Koordinaten */
    index* domelem;                  /** @brief Elemente, mit lokaler Knotennummerierung */
    index* c;                        /** @brief C Permutationsmatrix als Vektor */
    index* dombdry;                  /** @brief Knoten am Rand und Art der Randbed. */

    index nedgenodes;                /** @brief Anzahl Crosspoint+Edge Knoten -> Es gibt immer
                                         vier Crosspoints! Der domcoord Vektor ist ge-
                                         ordnet nach [crosspoints   edge nodes   interior nodes]' */

    index nfixed_loc;                /** @brief Anzahl der Gebietsrandpunkte in Domain */
    index* fixed_loc;                /** @brief Lokale Knotennummern der Knoten, welche auf
                                         dem Gebietsrand liegen */

    
    index neighbours[4];    /** @brief ranks of the neighbours [s e n w] if no neigh-
                                bour available -1 */
    index n_single_bdry[4]; /** @brief n_edgenodes (w.o. crosspoints) of each boundary 
                                    sortet [s e n w] */
    index n_cross_glob;     /** @brief nur auf rank 0, anzahl paarweise vers crosspoints */
    index* c_cross;         /** @brief nur auf rank 0, ordnet allen crosspoints von a
                                allen ranks die globale Knotennummer zu */
    index* cross_to_buf;    /** @brief nur auf rank 0, index fuer den buffer zur 
                                addition der crosspoints auf rank 0 */
    double* buf_cross;      /** @brief nur auf rank 0 buffer fuer gather crosspoints
                                Laenge des buffers ist 4*size */
    double* buf_cross_aggr; /** @brief nur auf rank 0 buffer um crosspoints 
                                aufzuaddieren, Laenge ist n_crosspoints_global */
                                
    double* buf_bdry[4];    /** @brief buffer um werte der boundaries empfangen zu 
                                koennen, 2_dimensionales array buf_brdy[0] fuer 
                                south bdry [1] fuer east usw. */
    bool black;             /** @brief meshes have a red black ordering, true is 
                                black, false is red */
    /* probably unneccessary
    index* map_south;
    index* map_east;
    index* map_north;
    index* map_west;
    */
} mesh_trans ;

mesh_trans* alloc_mesh_trans(index anz_dom, index dof);
mesh_trans* free_mesh_trans(mesh_trans* metra);

void mesh_buildRhs_loc(const mesh_trans *M, double *b, double (*fV)(double *, index),
                   double (*fN)(double *, index));
#endif
