#include <assert.h>
#include <time.h>

#include "hpc.h"

#ifndef GRID_M
#define GRID_M	4 
#endif
#ifndef GRID_N
#define GRID_N	3
#endif
#ifndef N_TESTRUNS_MAPPING
#define N_TESTRUNS_MAPPING 3
#endif


/**
  * @brief Test the mapping of vertices from locat meshes to the global mesh, by
  	   randomly picking local coordinates, determining the global index of
	   the vertex using the mapping and comparing the global with the local 
	   coordinates
  * @param globalMesh the global mesh
  * @param mapping Mapping object for the global mesh
  * @param procIndX x-index of the to be tested process in the MPI grid
  * @param procInd> y-index of the to be tested process in the MPI grid
  * @param rep number repetitions to be performed
  */
bool testVertexMapping(mesh *globalMesh, MeshMapping ***mapping, index procIndX,
		       index procIndY, int rep){
    bool testRes;
    index nLocalVert = mapping[procIndX][procIndY]->localMesh->ncoord;
    
    index randVal;
    double localCoords[2];
    double globalCoords[2];
    index mapIndGlbl;
    
    // perform specified number of repetitions
    for (int i=0; i<rep; ++i){
	// randomly select an index of a local vertex
	randVal = rand()%nLocalVert;
	// Get coords of randomly selected local vertex
	localCoords[0] = mapping[procIndX][procIndY]->localMesh->coord[randVal*2];
	localCoords[1] = mapping[procIndX][procIndY]->localMesh->coord[randVal*2+1];
	// Get global index through mapping
	mapIndGlbl = mapping[procIndX][procIndY]->vertexL2G[randVal];
	// Get global coordinates
	globalCoords[0] = globalMesh->coord[mapIndGlbl*2];
	globalCoords[1] = globalMesh->coord[mapIndGlbl*2+1];

	// determine test result
	testRes = (localCoords[0] == globalCoords[0] && localCoords[1] == globalCoords[1])?
		   1:0;

	printf("(%td,%td)|%-7td|(%1.3lf,%1.3lf)|(%1.3lf,%1.3lf)|%-4s\n", procIndX,
		procIndY, mapIndGlbl, localCoords[0], localCoords[1], globalCoords[0],
		globalCoords[1],(testRes)?"true":"false");

	if (!testRes){
	    break;
	}
    }

    return testRes;
}



/**
 * @brief compute the number of nodes for each local mesh
 * and write it to array
 *
 * @param buffer Array
 * @param globalMesh global Mesh
 * @param gridDimX x dimension of MPI Grid
 * @param gridDimY y dimension of MPI Grid
 */
void getTargetNodeCounts(index *buffer, mesh *globalMesh, index gridDimX,
                         index gridDimY) {
    // Compute the number of nodes for all local meshes
    for (index i = 0; i < gridDimY; ++i) {
        for (index j = 0; j < gridDimX; ++j) {
            buffer[i * gridDimX + j] =
                getSliceSize((index)sqrt((double)globalMesh->ncoord - 1),
                             gridDimY, i) +
                1;
            buffer[i * gridDimX + j] *=
                getSliceSize((index)sqrt((double)globalMesh->ncoord - 1),
                             gridDimX, j) +
                1;
        }
    }
}

/**
 * @brief Collect the numer of nodes of all local meshes and write the values to
 *	   annray
 * @param buffer array
 * @param mapping mesh mapping
 * @param gridDimX x dimension of MPI grid
 * @param gridDimY y dimension of MPI grid
 */
void collectLocalNodeCounts(index *buffer, MeshMapping ***mapping,
                            index gridDimX, index gridDimY) {
    for (int i = 0; i < gridDimY; ++i) {
        for (int j = 0; j < gridDimX; ++j) {
            buffer[i * gridDimX + j] = mapping[j][i]->localMesh->ncoord;
        }
    }
}

/**
 * @brief Compute the expected number of boundary elements for all lcoal meshes
 *	   and write the values to array
 * @param buffer array
 * @param mesh global Mesh
 * @param gridDimX x dimension of MPI grid
 * @param gridDimY y dimension of MPI grid
 */
void getTargetBndryCounts(index *buffer, mesh *mesh, index gridDimX,
                          index gridDimY) {
    index nEdgesX, nEdgesY;

    // First, compute expected number of edges for each local mesh
    for (index i = 0; i < gridDimY; ++i) {
        nEdgesY =
            getSliceSize((index)sqrt((double)mesh->ncoord - 1), gridDimY, i);

        for (index j = 0; j < gridDimX; ++j) {
            nEdgesX = getSliceSize((index)sqrt((double)mesh->ncoord - 1),
                                   gridDimX, j);

            /** determin expected number of boundary elements in dependence of
               the indices of the current local mesh.
             */
            buffer[i * gridDimX + j] = 0;

            if (i == 0 || i == gridDimY - 1) {
                buffer[i * gridDimX + j] = nEdgesX;
            }
            if (j == 0 || j == gridDimX - 1) {
                buffer[i * gridDimX + j] += nEdgesY;
            }
        }
    }
}

/**
  * @brief collect the number of boundary elements of all local meshes according
           according to the mapping and write the values to an array
  * @param buffer array
  * @param mapping mesh mapping
  * @param gridDimX x dimension of MPI grid
  * @param gridDimY y dimension of MPI grid
*/
void collectLocalBndryCounts(index *buffer, MeshMapping ***mapping,
                             index gridDimX, index gridDimY) {
    for (int i = 0; i < gridDimY; ++i) {
        for (int j = 0; j < gridDimX; ++j) {
            buffer[i * gridDimX + j] = mapping[j][i]->localMesh->nbdry;
        }
    }
}


/**
  * @brief Test for the mapping of vertices, elements and boundaries from the
	   local meshes to the global one. Tests are performed for the local
	   mesh of each process in the MPI grid
  * @param globalMesh the global unrefined mesh
  * @param gridDims dimensions of the MPI grid
*/
bool testMapping(mesh *globalMesh, index *gridDims){
    bool testResult = 0;
    
    // initial refinement of the mesh to make the split possible
    index initial_Refines = HPC_MAX(gridDims[0],gridDims[1]); 
    mesh *mesh_refined = mesh_initRefinement(globalMesh, initial_Refines);
    // spliet the mesh and thus create a mapping
    MeshMapping ***testMapping = mesh_split(mesh_refined,gridDims);
    
    printf("\n--- Test 2.1 Mapping of vertices ---\n"); 
    printf("%-5s|%-7s|%-13s|%-13s|%-4s\n", "proc", "map ind", "expect","actual","pass");
    printf("----------------------------------------------\n");
    // Test the mapping of the local mesh of each process in the MPI grid
    for (int i=0;i<gridDims[0]; ++i){
	for (int j=0; j<gridDims[1]; ++j){
	    testResult = testVertexMapping(mesh_refined, testMapping, i, j, N_TESTRUNS_MAPPING);
	    printf("----------------------------------------------\n");
	}
    }

    return testResult;
}



/**
  * @brief Test the partitioning of the mesh, cmoputed using the function
           mesh_split() in terms of the numer of nodes and boundary elements of
           each local mesh
  * @param globalMesh gloabl mesh
  * @param gridDims araay of dimensions of the MPI grid
  * @param numRefines number of grid refinements for wich the function shall be
            tested
*/
bool testNodePartitioning(mesh *globalMesh, index *gridDims, int numRefines) {	
    bool testResult = 1;
    
    mesh *mesh_current;
    mesh *mesh_previous;

    // initial refinement of the mesh to make the mapping possible according to
    // the mpi grid
    index initialRefines = HPC_MAX(gridDims[0], gridDims[1]);
    mesh_current = mesh_initRefinement(globalMesh, initialRefines);

    printf("\n%-4s|%-5s|%-14s|%-14s|%-8s|%-8s|%-5s\n", "ref", "(i,j)",
           "cnt nodes ref", "cnt nodes tst", "bdry ref", "bdry tst", "pass");
    printf("---------------------------------------------------------------\n");

    // Refine mesh iteratively and collect target and local node counts
    index refNodeCounts[gridDims[0] * gridDims[1]];
    index testNodeCounts[gridDims[0] * gridDims[1]];
    index refBndryCounts[gridDims[0] * gridDims[1]];
    index tstBndryCounts[gridDims[0] * gridDims[1]];
    MeshMapping ***testMapping;

    bool resultCurrent = 1;
    for (index i = initialRefines; i < numRefines; ++i) {
	// create mapping
        testMapping = mesh_split(mesh_current, gridDims);

        // Get raference node counts and node counts of local meshes
        getTargetNodeCounts(refNodeCounts, mesh_current, gridDims[0],
                            gridDims[1]);
        collectLocalNodeCounts(testNodeCounts, testMapping, gridDims[0],
                               gridDims[1]);

        // Get refenerence and local bndry edge counts
        getTargetBndryCounts(refBndryCounts, mesh_current, gridDims[0],
                             gridDims[1]);
        collectLocalBndryCounts(tstBndryCounts, testMapping, gridDims[0],
                                gridDims[1]);

        // Compare and print result
        index refCount, tstCount, refBndryCnts, tstBndryCnts;
        
	for (index k = 0; k < gridDims[1]; ++k) {
            for (index l = 0; l < gridDims[0]; ++l) {
                refCount = refNodeCounts[k * gridDims[0] + l];
                tstCount = testNodeCounts[k * gridDims[0] + l];
                refBndryCnts = refBndryCounts[k * gridDims[0] + l];
                tstBndryCnts = tstBndryCounts[k * gridDims[0] + l];
               
		resultCurrent = (refCount == tstCount &&
				refBndryCnts == tstBndryCnts)?1:0;		

		printf("%-4td|(%td,%td)|%-14td|%-14td|%-8td|%-8td|%-5s\n", i, l,
                       k, refCount, tstCount, refBndryCnts, tstBndryCnts,
		       resultCurrent? "true":"FALSE");
		
		// Set test Result to false if current test did not pass
		if (resultCurrent == 0){
		    testResult = 0;
		}
            }
        }
        printf(
            "---------------------------------------------------------------"
            "\n");
        // next refinement
        mesh_previous = mesh_current;
        mesh_current = mesh_refine(mesh_previous);
        mesh_getEdge2no(mesh_current->nelem, mesh_current->elem,
                        &mesh_current->nedges, &mesh_current->edge2no);
        mesh_current->fixed =
            mesh_getFixed(mesh_current->ncoord, mesh_current->bdry,
                          mesh_current->nbdry, &mesh_current->nfixed);

        // delete old mapping
        delete2DMeshMapping(testMapping, gridDims[0]);

        // free old mesh
        mesh_free(mesh_previous);
    }

    // free final mesh
    mesh_free(mesh_current);

    return testResult;
}

int main() {
    srand(time(NULL));

    // Init global result vector to be negative 
    bool resultTotal = 1;
    bool resultCurrent;

    printf("\n=== Start test_mesh_split ===\n");

    char fname[32] = "../Problem/problem1";
    mesh *m1 = mesh_load(fname);
    
    mesh_getEdge2no(m1->nelem, m1->elem, &m1->nedges, &m1->edge2no);
    m1->fixed = mesh_getFixed(m1->ncoord, m1->bdry, m1->nbdry, &m1->nfixed);
  
    printf("\n--- Test 1: Partitioning of Nodes ---\n");
    index dims[2] = {GRID_N,GRID_M};

    resultCurrent = testNodePartitioning(m1, dims, 8);

    // rewrite total test result only if the test didn't already fail 
    if (resultTotal) resultTotal = resultCurrent;

    printf("\n--- Test 2: Mapping ---\n");
    resultCurrent = testMapping(m1, dims);

    if (resultTotal) resultTotal = resultCurrent; 

    printf("\n--- Total Result: %-4s---\n",resultCurrent?"PASS":"FAIL");

    printf("=== End test_mesh_split ===\n");
    mesh_free(m1);
}
