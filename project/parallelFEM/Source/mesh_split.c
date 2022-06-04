#include <assert.h>

#include "hpc.h"

/**
 * @brief Get x/y-coordinate of vertex i.
 *
 * @param m Mesh
 * @param i Index
 * @param j 0 for x-, 1 for y-coordinate
 * @return double x-coordinate
 */
double getCoord(mesh *m, index i, index j) {
    assert(j == 0 || j == 1);
    return m->coord[2 * i + j];
}

/**
 * @brief Get the x/y-index of the x/y-coordinate.
 *
 * @param gMeshDim Number of vertices in x/y-direction
 * @param pxy x/y-coordinate
 * @return index x/y-index
 */
index getIndex(index gMeshDim, double pxy) {
    // Compute the size of the domain each process gets
    double h = 1. / (gMeshDim - 1);
    index xy = (index)(pxy / h);
    assert(xy >= 0);
    assert(xy < gMeshDim);
    return xy;
}

/**
 * @brief Returns local index of global vertex with globalIndex.
 *
 * @param mapping
 * @param globalIndex
 * @return index
 */
index vertexG2L(MeshMapping *mapping, index globalIndex) {
    for (index i = 0; i < mapping->localMesh->ncoord; i++) {
        if (globalIndex == mapping->vertexL2G[i]) {
            return i;
        }
    }
    printf("Local vertex index for global vertex index %zu not found!\n", globalIndex);
    abort();
}

/**
 * @brief Inserts coordinates of global vertex i into mapping->localMesh
 * coordinates as local vertex j. Saves the local to global vertex mapping j->i.
 *
 * @param globalMesh
 * @param i
 * @param mapping
 * @param j
 */
void insertCoord(mesh *globalMesh, index i, MeshMapping *mapping, index j) {
    assert(i >= 0);
    assert(i < globalMesh->ncoord);
    assert(j >= 0);
    assert(j < mapping->localMesh->ncoord);
    mapping->localMesh->coord[2 * j] = globalMesh->coord[2 * i];
    mapping->localMesh->coord[2 * j + 1] = globalMesh->coord[2 * i + 1];
    mapping->vertexL2G[j] = i;
}

/**
 * @brief Inserts global elem (triangle) i into mapping->localMesh as elem
 * (triangle) j.
 *
 * @param globalMesh
 * @param i
 * @param mapping
 * @param j
 */
void insertElem(mesh *globalMesh, index i, MeshMapping *mapping, index j) {
    assert(globalMesh != NULL);
    assert(i >= 0);
    assert(i < globalMesh->nelem);
    assert(mapping->localMesh != NULL);
    assert(j >= 0);
    assert(j < mapping->localMesh->nelem);
    for (index k = 0; k < 7; k++) {
        if (k < 3) {
            mapping->localMesh->elem[7 * j + k] =
                vertexG2L(mapping, globalMesh->elem[7 * i + k]);
        } else {
            mapping->localMesh->elem[7 * j + k] = globalMesh->elem[7 * i + k];
        }
    }
    mapping->elemL2G[j] = i;
}

/**
 * @brief Inserts global bdry edge i into mapping->localMesh as bdry edge j.
 *
 * @param globalMesh
 * @param i
 * @param mapping
 * @param j
 */
void insertBdry(mesh *globalMesh, index i, MeshMapping *mapping, index j) {
    assert(globalMesh != NULL);
    assert(i >= 0);
    assert(i < globalMesh->nbdry);
    assert(mapping->localMesh != NULL);
    assert(j >= 0);
    assert(j < mapping->localMesh->nbdry);
    for (index k = 0; k < 4; k++) {
        if (k < 2) {
            mapping->localMesh->bdry[4 * j + k] =
                vertexG2L(mapping, globalMesh->bdry[4 * i + k]);
        } else {
            mapping->localMesh->bdry[4 * j + k] = globalMesh->bdry[4 * i + k];
        }
    }
    mapping->bdryL2G[j] = i;
}

/**
 * @brief Get x/y-index of x/y-coordinate.
 *
 * @param gMeshDim
 * @param gridDimXY
 * @param pxy
 * @return index
 */
index getGridIndexOfCoord(index gMeshDim, index gridDimXY, double pxy) {
    index xy = getIndex(gMeshDim, pxy);
    return getSliceIndex(gMeshDim - 1, gridDimXY, xy);
}

/**
 * @brief Get grid index of vertex
 *
 * @param globalMesh
 * @param gMeshDim
 * @param gridDimXY
 * @param i Index of vertex
 * @param j 0 for x-, 1 for y-index
 * @return index
 */
index getGridIndexOfVertex(mesh *globalMesh,                 //
                           index gMeshDim, index gridDimXY,  //
                           index i, index j) {
    assert(j == 0 || j == 1);
    double pxy = getCoord(globalMesh, i, j);
    return getGridIndexOfCoord(gMeshDim, gridDimXY, pxy);
}

/**
 * @brief Get grid index of elem vertex
 *
 * @param globalMesh
 * @param gMeshDim
 * @param gridDimX
 * @param i Index of elem
 * @param j Index of elem vertex; Must be 0, 1 or 2
 * @param k 0 for x-, 1 for y-index
 * @return index
 */
index getGridIndexOfElemVertex(mesh *globalMesh,                 //
                               index gMeshDim, index gridDimXY,  //
                               index i, index j, index k) {
    assert(j == 0 || j == 1 || j == 2);
    assert(k == 0 || k == 1);
    double pxy = getCoord(globalMesh, globalMesh->elem[7 * i + j], k);
    return getGridIndexOfCoord(gMeshDim, gridDimXY, pxy);
}

/**
 * @brief Get grid index of bdry vertex
 *
 * @param globalMesh
 * @param gMeshDim
 * @param gridDimXY
 * @param i Index of bdry edge
 * @param j Index of bdry edge vertex; Must be 0 or 1
 * @param k 0 for x-, 1 for y-index
 * @return index
 */
index getGridIndexOfBdryVertex(mesh *globalMesh,                 //
                               index gMeshDim, index gridDimXY,  //
                               index i, index j, index k) {
    assert(j == 0 || j == 1);
    assert(k == 0 || k == 1);
    double pxy = getCoord(globalMesh, globalMesh->bdry[4 * i + j], k);
    return getGridIndexOfCoord(gMeshDim, gridDimXY, pxy);
}

/**
 * @brief This function returns the index local mesh index of an E vertex.
 *
 * @param xy The global x/y-index of the vertex
 * @param offsetKL Offset of the slice to which vertex with index x/y belongs to
 * @param offsetDir Offset depending on the cardinal direction
 * @return index
 */
index getLocalBorderIndex(index xy, index offsetKL, index offsetDir) {
    index i = 4 + (xy - 1) - offsetKL + offsetDir;
    // printf("%s: %zd\n", lMeshDimXY == 0 ? "s/w" : "n/e", i);
    return i;
}

/**
 * @brief Splits the mesh for the given grid dimension
 *
 * @param globalMesh The global mesh
 * @param gridDims The dimensions of the grid (MPI_Cart_get(..., gridDims, ...))
 * @return MeshMapping *** This is a 2D-array of dimensions gridDims
 */
MeshMapping ***mesh_split(mesh *globalMesh, int gridDims[2]) {
    // Dimensions of the 2d MPI grid
    index gridDimX = gridDims[0];
    index gridDimY = gridDims[1];
    MeshMapping ***mapping = new2DMeshMapping(gridDimX, gridDimY);

    // Global mesh variables
    index gncoord = globalMesh->ncoord;
    index gnelem = globalMesh->nelem;
    index gnbdry = globalMesh->nbdry;
    double *gcoord = globalMesh->coord;
    index gMeshDim = (index)sqrt((double)gncoord);
    assert(gMeshDim * gMeshDim == gncoord);

    // Allocate memory for indices, necessary for ordering
    index indicesV[gridDimX][gridDimY];
    index indicesI[gridDimX][gridDimY];
    index indicesE[gridDimX][gridDimY];

    // Compute local mesh variables and allocate mesh and mapping memory
    for (index k = 0; k < gridDimX; k++) {
        for (index l = 0; l < gridDimY; l++) {
            // Compute local mesh dimensions
            // Plus/minus 1 to account for vertices at start and end
            index lMeshDimX = getSliceSize(gMeshDim - 1, gridDimX, k) + 1;
            index lMeshDimY = getSliceSize(gMeshDim - 1, gridDimY, l) + 1;

	    DEBUG_PRINT("[%td,%td] lMeshDimX=%td, lMeshDimY=%td\n",k,l,
			lMeshDimX, lMeshDimY);

            // Compute local ncoord, nelem, nbdry
            index lncoord = lMeshDimX * lMeshDimY;
            index lnelem = 2 * (lMeshDimX - 1) * (lMeshDimY - 1);
            index lnedges = lMeshDimX * (lMeshDimY - 1)           // | edges
                            + lMeshDimY * (lMeshDimX - 1)         // - edges
                            + (lMeshDimX - 1) * (lMeshDimY - 1);  // / edges
            index lnbdry = 0;
            if (k == 0) {
                lnbdry += lMeshDimY - 1;
            }
	    if (k +1 == gridDimX){
		lnbdry += lMeshDimY -1;
	    }
            if (l == 0) {
                lnbdry += lMeshDimX - 1;
            }
	    if (l + 1 == gridDimY){
		lnbdry += lMeshDimX -1;
	    }

            // Allocate mesh and mapping memory
            mesh *localMesh = mesh_alloc(lncoord, lnelem, lnbdry);
            localMesh->nedges = lnedges;
            mapping[k][l] = newMeshMapping(localMesh, gncoord, lMeshDimX, lMeshDimY);

            // The V four nodes are first,
            // after that we have the 2 * (lMeshDimX - 2) + 2 * (lMeshDimY - 2)
            // E nodes, last come the remaining I nodes.
            indicesV[k][l] = 0;
            indicesE[k][l] = 4;
            indicesI[k][l] = 2 * (lMeshDimX - 2) + 2 * (lMeshDimY - 2) + 4;
        }
    }

    // Assign vertices from global mesh to local meshes
    for (index i = 0; i < gncoord; i++) {
        // (x,y) are the plane indices of vertex i
        index x = getIndex(gMeshDim, getCoord(globalMesh, i, 0));
        index y = getIndex(gMeshDim, getCoord(globalMesh, i, 1));

        // (k,l) are grid indices vertex i=(x,y) belongs to
        index k = getGridIndexOfVertex(globalMesh, gMeshDim, gridDimX, i, 0);
        index l = getGridIndexOfVertex(globalMesh, gMeshDim, gridDimY, i, 1);

        /**
         * If (x == offsetK and y == offsetL), then the current point is a V
         * vertex.
         * If the above does not hold, but (x == offsetK or y == offsetL) is
         * true, then the current point is an E vertex.
         * Otherwise it is an I vertex.
         */
        index offsetK = getSliceOffset(gMeshDim - 1, gridDimX, k);
        index offsetL = getSliceOffset(gMeshDim - 1, gridDimY, l);
        index k1 = k - 1;
        index l1 = l - 1;
        if (x == offsetK && y == offsetL) {  // V
            // Numbering:
            // 23
            // 01
            // We may have to insert the vertex into 4 meshes,
            // i.e. in all directions
            // north-east
            if (k < gridDimX && l < gridDimY) {
                indicesV[k][l]++;
                insertCoord(globalMesh, i, mapping[k][l], 0);
            }
            // north-west
            if (k1 >= 0 && l < gridDimY) {
                indicesV[k1][l]++;
                insertCoord(globalMesh, i, mapping[k1][l], 1);
            }
            // south-east
            if (k < gridDimX && l1 >= 0 && y == offsetL) {
                indicesV[k][l1]++;
                insertCoord(globalMesh, i, mapping[k][l1], 2);
            }
            // south-wests
            if (k1 >= 0 && l1 >= 0 && x == offsetK && y == offsetL) {
                indicesV[k1][l1]++;
                insertCoord(globalMesh, i, mapping[k1][l1], 3);
            }
        } else if (x == offsetK) {  // E, on a | border
            // We may have to insert the vertex into 2 meshes,
            // i.e. east and west
            if (k < gridDimX && l < gridDimY) {
                indicesE[k][l]++;
                index offsetDir = 2 * (mapping[k][l]->lMeshDimX - 2);
                insertCoord(globalMesh, i, mapping[k][l], getLocalBorderIndex(y, offsetL, offsetDir));
            }
            if (k1 >= 0 && l < gridDimY) {
                indicesE[k1][l]++;
                index offsetDir = 2 * (mapping[k1][l]->lMeshDimX - 2) + mapping[k1][l]->lMeshDimY - 2;
                insertCoord(globalMesh, i, mapping[k1][l], getLocalBorderIndex(y, offsetL, offsetDir));
            }
        } else if (y == offsetL) {  // E, on a - border
            // We may have to insert the vertex into 2 meshes,
            // i.e. north and south
            index offsetDir = 0;
            if (k < gridDimX && l < gridDimY) {
                indicesE[k][l]++;
                insertCoord(globalMesh, i, mapping[k][l], getLocalBorderIndex(x, offsetK, offsetDir));
            }
            if (k < gridDimX && l1 >= 0) {
                indicesE[k][l1]++;
                offsetDir += mapping[k][l1]->lMeshDimX - 2;
                insertCoord(globalMesh, i, mapping[k][l1], getLocalBorderIndex(x, offsetK, offsetDir));
            }
        } else {  // I
            assert(k < gridDimX);
            assert(l < gridDimY);
            insertCoord(globalMesh, i, mapping[k][l], indicesI[k][l]++);
        }
    }

    // Asserting that we got all elements
    for (index k = 0; k < gridDimX; k++) {
        for (index l = 0; l < gridDimY; l++) {
            assert(indicesV[k][l] == 4);
            assert(indicesE[k][l] == 2 * (mapping[k][l]->lMeshDimX - 2) + 2 * (mapping[k][l]->lMeshDimY - 2) + 4);
            assert(indicesI[k][l] == mapping[k][l]->localMesh->ncoord);
        }
    }

    // Assign elements to local meshes
    index indices[gridDimX][gridDimY];
    memset(indices, 0, gridDimX * gridDimY * sizeof(index));
    for (index i = 0; i < gnelem; i++) {
        index k1 = getGridIndexOfElemVertex(globalMesh, gMeshDim, gridDimX, i, 0, 0);
        index k2 = getGridIndexOfElemVertex(globalMesh, gMeshDim, gridDimX, i, 1, 0);
        index k3 = getGridIndexOfElemVertex(globalMesh, gMeshDim, gridDimX, i, 2, 0);
        index k = HPC_MIN(k1, HPC_MIN(k2, k3));

        index l1 = getGridIndexOfElemVertex(globalMesh, gMeshDim, gridDimY, i, 0, 1);
        index l2 = getGridIndexOfElemVertex(globalMesh, gMeshDim, gridDimY, i, 1, 1);
        index l3 = getGridIndexOfElemVertex(globalMesh, gMeshDim, gridDimY, i, 2, 1);
        index l = HPC_MIN(l1, HPC_MIN(l2, l3));

        insertElem(globalMesh, i, mapping[k][l], indices[k][l]++);
    }

    // Asserting that we got all elements
    for (index k = 0; k < gridDimX; k++) {
        for (index l = 0; l < gridDimY; l++) {
            assert(indices[k][l] == mapping[k][l]->localMesh->nelem);
        }
    }

    // Assign boundary edges to local meshes
    memset(indices, 0, gridDimX * gridDimY * sizeof(index));
    for (index i = 0; i < gridDimX; i++){
        for (index j = 0; j < gridDimY; j++){
            DEBUG_PRINT("mapping[%zd][%zd]->localMesh->nbdry=%zd\n", i, j,
                        mapping[i][j]->localMesh->nbdry);
	 }
    }

    for (index i = 0; i < gnbdry; i++) {
        index k1 = getGridIndexOfBdryVertex(globalMesh, gMeshDim, gridDimX, i, 0, 0);
        index l1 = getGridIndexOfBdryVertex(globalMesh, gMeshDim, gridDimY, i, 0, 1);
        index k2 = getGridIndexOfBdryVertex(globalMesh, gMeshDim, gridDimX, i, 1, 0);
        index l2 = getGridIndexOfBdryVertex(globalMesh, gMeshDim, gridDimY, i, 1, 1);
        DEBUG_PRINT("%3zd | (%zd,%zd)-(%zd,%zd) | ", i, k1, l1, k2, l2);
        DEBUG_PRINT("(%5.3lf,%5.3lf)-(%5.3lf,%5.3lf) | ",
                    getCoord(globalMesh, globalMesh->bdry[4 * i], 0),
                    getCoord(globalMesh, globalMesh->bdry[4 * i], 1),
                    getCoord(globalMesh, globalMesh->bdry[4 * i + 1], 0),
                    getCoord(globalMesh, globalMesh->bdry[4 * i + 1], 1));

        index k = HPC_MIN(k1, k2);
        index l = HPC_MIN(l1, l2);

        /**
         *      north
         *      ^__
         *      |  |
         * west |__| > east
         *      south
         */
        if (k < gridDimX && l < gridDimY) {  // west or south
            DEBUG_PRINT("%3s | indices[%zd][%zd]=%2zd | nbdry=%zd\n", "w/s", k,
                        l, indices[k][l], mapping[k][l]->localMesh->nbdry);
            insertBdry(globalMesh, i,  //
                       mapping[k][l], indices[k][l]++);
        } else if (l == gridDimY) {  // north
            DEBUG_PRINT("%3s | indices[%zd][%zd]=%2zd | nbdry=%zd\n", "n", k,
                        l - 1, indices[k][l - 1],
                        mapping[k][l - 1]->localMesh->nbdry);
            insertBdry(globalMesh, i,  //
                       mapping[k][l - 1], indices[k][l - 1]++);
        } else if (k == gridDimX) {  // east
            DEBUG_PRINT("%3s | indices[%zd][%zd]=%2zd | nbdry=%zd\n", "e",
                        k - 1, l, indices[k - 1][l],
                        mapping[k - 1][l]->localMesh->nbdry);
            insertBdry(globalMesh, i,  //
                       mapping[k - 1][l], indices[k - 1][l]++);
        } else {
            printf("Something went wrong! :(\n");
            abort();
        }
    }

    // Asserting that we got all boundaries
    for (index k = 0; k < gridDimX; k++) {
        for (index l = 0; l < gridDimY; l++) {
            if (indices[k][l] != mapping[k][l]->localMesh->nbdry) {
                DEBUG_PRINT("indices[%zd][%zd]=%zd, ", k, l, indices[k][l]);
                DEBUG_PRINT("mapping[%zd][%zd]->localMesh->nbdry=%zd\n", k, l,
                            mapping[k][l]->localMesh->nbdry);
            }
            assert(indices[k][l] == mapping[k][l]->localMesh->nbdry);
        }
    }
    DEBUG_PRINT("%s\n", "Passt");

    return mapping;
}
