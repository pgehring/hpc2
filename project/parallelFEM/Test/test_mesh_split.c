#include <assert.h>
#include "hpc.h"
#include <math.h>

#ifndef GRID_M
#define GRID_M 4
#endif
#ifndef GRID_N
#define GRID_N 3
#endif
#ifndef N_TESTRUNS_MAPPING
#define N_TESTRUNS_MAPPING 3
#endif

/**
  * @brief function to determine the expected number of nodes  in given
	   local mesh
  * @param globalMesh the global Mesh
  * @param gridIndX x-Index of the process in the MPI grid
  * @param gridIndY y-Index of the process in the MPI grid
  * @param gridDims dimensions of the two dimensional MPI grid
  * @return expected number of nodes  in the given local mesh
*/
index getTargetNodeCount(mesh *globalMesh, int gridIndX, int gridIndY,
			 int *gridDims){
    index localXDim = getSliceSize((index)sqrt((double)globalMesh->ncoord -1),
				   (index)gridDims[0], gridIndX) + 1;
    index localYDim = getSliceSize((index)sqrt((double)globalMesh->ncoord - 1),
				   (index)gridDims[1], gridIndY) + 1;
    return localXDim * localYDim;
}


/**
  * @breif wrapper function that yields the number of nodes in the local mesh
  * @param mapping pointer to two dimensional array of mapping structs, i.e. the
		   global mapping
  * @param gridIndX x-Index in the MPI grid
  * @param gridIndY y-Index in the MPI grid
  * @param gridDims dimensions of the two dimensional MPI grid
  * @return number of nodes in the local mesh
*/
index getLocalNodeCount(MeshMapping ***mapping, int gridIndX, int gridIndY,
			    int *gridDims){
    return mapping[gridIndX][gridIndY]->localMesh->ncoord;
}


/**
  * @brief function to determine the expected number of elements  in given
	   local mesh
  * @param globalMesh the global Mesh
  * @param gridIndX x-Index of the process in the MPI grid
  * @param gridIndY y-Index of the process in the MPI grid
  * @param gridDims dimensions of the two dimensional MPI grid
  * @return expected number of elements in the given local mesh
*/
index getTargetElementCount(mesh *globalMesh, int gridIndX, int gridIndY,
			    int *gridDims){
    index glblMeshDim = (index)sqrt((double)globalMesh->ncoord);
    index nEdgesX = getSliceSize(glblMeshDim-1, gridDims[0], gridIndX);
    index nEdgesY = getSliceSize(glblMeshDim-1, gridDims[1], gridIndY);
    
    return 2*nEdgesX*nEdgesY;
}


/**
  * @breif wrapper function that yields the number of elements in the local mesh
  * @param mapping pointer to two dimensional array of mapping structs, i.e. the
		   global mapping
  * @param gridIndX x-Index in the MPI grid
  * @param gridIndY y-Index in the MPI grid
  * @param gridDims dimensions of the two dimensional MPI grid
  * @return number of elements in the local mesh
*/
index getLocalElemCount(MeshMapping ***mapping, int gridIndX, int gridIndY,
			int *gridDims){
    return mapping[gridIndX][gridIndY]->localMesh->nelem;
}

/**
  * @brief function to determine the expected number of boundarie nodes in given
	   local mesh
  * @param globalMesh the global Mesh
  * @param gridIndX x-Index of the process in the MPI grid
  * @param gridIndY y-Index of the process in the MPI grid
  * @param gridDims dimensions of the two dimensional MPI grid
  * @return expected number of boundary nodes in the given local mesh
*/
index getTargetBndryCount(mesh *globalMesh, int gridIndX, int gridIndY,
			  int *gridDims){
    index targetBdryCnt = 0;
    
    index nEdgesX = getSliceSize((index)sqrt((double)globalMesh->ncoord - 1),
				 gridDims[0], gridIndX);
    index nEdgesY = getSliceSize((index)sqrt((double)globalMesh->ncoord - 1),
				 gridDims[1], gridIndY);
   

    if (gridIndX == 0 || gridIndX == gridDims[0] -1){
	targetBdryCnt = nEdgesY;
    }
    if (gridIndY == 0 || gridIndY == gridDims[1] -1){
	targetBdryCnt += nEdgesX;
    }
    return targetBdryCnt;
}

/**
  * @breif wrapper function that yields the number of boundary nodes in the local mesh
  * @param mapping pointer to two dimensional array of mapping structs, i.e. the
		   global mapping
  * @param gridIndX x-Index in the MPI grid
  * @param gridIndY y-Index in the MPI grid
  * @param gridDims dimensions of the two dimensional MPI grid
  * @return number of boundary nodes in the local mesh
*/
index getLocalBndryCount(MeshMapping ***mapping, int gridIndX, int gridIndY,
			 int *gridDims){
    return mapping[gridIndX][gridIndY]->localMesh->nbdry;
}

/**
  * @brief function for testing the partitioning of nodes performed by splitting
	   the global mesh using function mesh_split. The test is done with
	   respect to the number node in the local meshes, the number elements,
	   and boundary nodes.
    @param globalMesh the global mesh struct
    @param gridDims dimensions of the two dimensional MPI grid
    @param numRefines number of refinemenets of the initial mesh that shall be
	   performed during the test
    @return testResult boolean value representing the result of the test
*/
bool testNodePartitioning(mesh *globalMesh, int gridDims[2], int numRefines){
    bool testResult = 1;

    mesh *mesh_current;
    mesh *mesh_previous;

    // initial refinement of the mesh to make the mapping possible according to
    // the mpi grid
    index initialRefines = HPC_MAX(gridDims[0], gridDims[1]);
    mesh_current = mesh_initRefinement(globalMesh, initialRefines);

    printf("\n%-4s|%-5s|%-14s|%-14s|%-8s|%-8s|%-8s|%-8s|%-5s\n", "ref", "(x,y)",
           "cnt nodes ref", "cnt nodes tst","elem ref", "elem tst", "bdry ref",
	   "bdry tst", "pass");
    printf("----------------------------------------------------------------------------------\n");
    
    index refNodeCnt, tstNodeCnt, refBndryCnt, tstBndryCnt, refElemCnt, tstElemCnt;
    MeshMapping ***testMapping;
    bool resultCurrent = 1;

    // iteratively refine the mesh and do the tests
    for (index i = initialRefines; i < numRefines; ++i) {
	printf("global Meshdim: %td\n",(index)sqrt((double)mesh_current->ncoord));

	// create mapping
        testMapping = mesh_split(mesh_current, gridDims);

	// go through the MPI grid column wise
	for (index k = 0; k < gridDims[1]; ++k) {

	    for (index l = 0; l < gridDims[0]; ++l) {
		// get reference and local node count for current process
		refNodeCnt = getTargetNodeCount(mesh_current, l, k, gridDims);
		tstNodeCnt = getLocalNodeCount(testMapping, l, k, gridDims);

		// get reference and local element count for current process
		refElemCnt = getTargetElementCount(mesh_current, l, k, gridDims);
		tstElemCnt =  getLocalElemCount(testMapping, l, k, gridDims);

		// get reference and local bndry count for current process
		refBndryCnt = getTargetBndryCount(mesh_current, l, k, gridDims);
		tstBndryCnt = getLocalBndryCount(testMapping, l, k, gridDims);

		// determine intermediary test result
		resultCurrent = (refNodeCnt == tstNodeCnt &&
				 refBndryCnt == tstBndryCnt &&
				 refElemCnt == tstElemCnt)?1:0;


		printf("%-4td|(%td,%td)|%-14td|%-14td|%-8td|%-8td|%-8td|%-8td|%-5s\n",
		       i, l, k, refNodeCnt, tstNodeCnt, refElemCnt, tstElemCnt, 
		       refBndryCnt, tstBndryCnt,
                       resultCurrent ? "true" : "FALSE");

		// Set test Result to false if current test did not pass
                if (resultCurrent == 0) {
                    testResult = 0;
                }
	    }
	}

        printf(
            "----------------------------------------------------------------------------------"
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

    return testResult;
}

/**
  * @brief compute root mean square error of the entries of two arrays of doubles
  * @param refEntries array of reference entries
  * @param tstEntries array of test entries
  * @param number of entries in the arrays (for both the same)
*/
double computeRMSE(double *refEntries, double *tstEntries, int nValues){
    double rmse = 0;
    double diff; 

    for (int i=0; i<nValues; ++i){
	diff = refEntries[i] - tstEntries[i];
	rmse += pow(diff,2.0);
    }

    rmse = sqrt(rmse);

    return rmse;
}


/**
  * @brief Get the coordinates of nodes of a boundary edge by a boundary index in 
	   the local scope of the process
  * @param buffer memory area to store the coordinate values with length 4
	   [(x1,y1),(x2,y2)]
  * @param mapping the global mapping (pointer to two dimensional array of mapping
	   structs)
  * @param bndryInd index of the boundary edge in the local scope of the process
  * @param procIndX x-coordinate of the process in the MPI Grid
  * @param procIndY y-coordinate of the process in the MPI Grid
*/
void getLocalNodeCoordsForBoundary(double *buffer, MeshMapping ***mapping,
				   index bndryInd, index procIndX, index procIndY){
    index node0, node1;
    // get local node indices
    node0 = mapping[procIndX][procIndY]->localMesh->bdry[bndryInd*4];
    node1 = mapping[procIndX][procIndY]->localMesh->bdry[bndryInd*4+1];

    // insert coordinates of the first node
    buffer[0] = mapping[procIndX][procIndY]->localMesh->coord[node0*2];
    buffer[1] = mapping[procIndX][procIndY]->localMesh->coord[node0*2+1];
    // insert coordinates of the second node
    buffer[2] = mapping[procIndX][procIndY]->localMesh->coord[node1*2];
    buffer[3] = mapping[procIndX][procIndY]->localMesh->coord[node1*2+1];
}

/**
  * @brief Get the coordinates of nodes of a boundary edge by its global index
  * @param buffer memory area to store the coordinate values with length 4
	   [(x1,y1),(x2,y2)]
  * @param globalMesh pointer to the global mesh struct
  * @param bndryInd global boundary index
  * @procIndX x-coordinate of the process in the MPI Grid
  * @procIndY y-coordinate of the process in the MPI Grid
*/
void getGlobalNodeCoordsForBoundary(double *buffer, mesh *globalMesh,
				   index bndryInd, index procIndX, index procIndY){
    index node0, node1;
    // get local node indices
    node0 = globalMesh->bdry[bndryInd*4];
    node1 = globalMesh->bdry[bndryInd*4+1];

    // insert coordinates of first node
    buffer[0] = globalMesh->coord[node0*2];
    buffer[1] = globalMesh->coord[node0*2+1];
    // insert coordinates of second node
    buffer[2] = globalMesh->coord[node1*2];
    buffer[3] = globalMesh->coord[node1*2+1];
}


/**
  * @brief function to test the mapping of the boundary edges from the local to
	   the global scope. The test is performed by randomly selecting a local
	   boundary index, getting the global index using the mapping and comapring
	   the coordinates of the corresponding nodes of the boundary edge.
  * @param globalMesh the global mesh struct
  * @param mapping pointer to two dimensional array of mapping structs i.e. the
	   global mapping
  * @param procIndX x-coordinate of the process
  * @param procIndY y-coordinate of the process
  * @param rep number of repetitions to be performed
  * @return root mean square error of the compared node coordinates
*/
double testBoundaryMapping(mesh *globalMesh, MeshMapping ***mapping, index procIndX,
			   index procIndY, int rep){
    double rmseTotal = 0;
    index nLocalBndry = mapping[procIndX][procIndY]->localMesh->nbdry;

    // If the current local mesh does not include boundary edges, return -1
    if (nLocalBndry == 0){
	return -1.0;
    }

    index randVal, glblBndryInd;
    double localCoords[4], globalCoords[4], rmseCurrent;

    // perform specified number of repetitions
    for (int i=0; i<rep; ++i){
	//randomly select an index of a local vertex
	randVal = rand()%nLocalBndry;

	// get local boundary coords for the boundary edge at local ind randVal
	getLocalNodeCoordsForBoundary(localCoords, mapping, randVal, procIndX, procIndY);
	
	// get global boundary index
	glblBndryInd = mapping[procIndX][procIndY]->bdryL2G[randVal];
	getGlobalNodeCoordsForBoundary(globalCoords, globalMesh, glblBndryInd,
				       procIndX, procIndY);

	// Compute RMSE
	rmseCurrent = computeRMSE(globalCoords, localCoords, 4);
	rmseTotal = sqrt(pow(rmseTotal,2.0)+pow(rmseCurrent,2.0));
	
	DEBUG_PRINT("\nDEBUG BOUNDARY MAPPING\n","");
	DEBUG_PRINT("ref: [(%lf,%lf),(%lf,%lf)], tst: [(%lf,%lf),(%lf,%lf)] rmse: %lf\n",
		    globalCoords[0], globalCoords[1], globalCoords[2], globalCoords[3],
		    localCoords[0], localCoords[1], localCoords[2],localCoords[3],
		    rmseCurrent);
    }
    return rmseTotal;
}


/**
  * @brief Get the coordinates of nodes of an element by an element index in 
	   the local scope of the process
  * @param buffer memory area to store the coordinate values with length 6
	   [(x1,y1),(x2,y2),(x3,y3)]
  * @param mapping the global mapping (pointer to two dimensional array of mapping
	   structs)
  * @param elemInd index of the element in the local scope of the process
  * @param procIndX x-coordinate of the process in the MPI Grid
  * @param procIndY y-coordinate of the process in the MPI Grid
*/
void getLocalNodeCoordsForElement(double *buffer, MeshMapping ***mapping,
				  index elemInd, index procIndX, index procIndY){
    index localNodeInds[3];

    // insert local node indices
    localNodeInds[0] = mapping[procIndX][procIndY]->localMesh->elem[elemInd*7];
    localNodeInds[1] = mapping[procIndX][procIndY]->localMesh->elem[elemInd*7+1];
    localNodeInds[2] = mapping[procIndX][procIndY]->localMesh->elem[elemInd*7+2];
    
    // insert coordinates matching the local node indices
    for (int i=0; i<3; ++i){
	buffer[i*2] = mapping[procIndX][procIndY]->localMesh->coord[localNodeInds[i]*2];
	buffer[i*2+1] = mapping[procIndX][procIndY]->localMesh->coord[localNodeInds[i]*2+1];
    }
}


/**
  * @brief Get the coordinates of nodes of an element by its global index
  * @param buffer memory area to store the coordinate values with length 6
	   [(x1,y1),(x2,y2),(x3,y3)]
  * @param globalMesh pointer to the global mesh struct
  * @param glblElemInd global element index
  * @procIndX x-coordinate of the process in the MPI Grid
  * @procIndY y-coordinate of the process in the MPI Grid
*/
void getGlobalNodeCoordsForElement(double *buffer, mesh *globalMesh, 
				  index glblElemInd, index procIndX, index procIndY){
    index globalNodeInds[3];
    
    // inser global node indices
    globalNodeInds[0] = globalMesh->elem[glblElemInd*7];
    globalNodeInds[1] = globalMesh->elem[glblElemInd*7+1];
    globalNodeInds[2] = globalMesh->elem[glblElemInd*7+2];

    // inser coordinates matching the global node indices
    for (int i=0; i<3; ++i){
	buffer[i*2] = globalMesh->coord[globalNodeInds[i]*2];
	buffer[i*2+1] = globalMesh->coord[globalNodeInds[i]*2+1];
    }
}


/**
  * @brief function to test the mapping of the elements from the local to
	   the global scope. The test is performed by randomly selecting a local
	   element index, getting the global index using the mapping and comapring
	   the coordinates of the corresponding nodes of the element.
  * @param globalMesh the global mesh struct
  * @param mapping pointer to two dimensional array of mapping structs i.e. the
	   global mapping
  * @param procIndX x-coordinate of the process
  * @param procIndY y-coordinate of the process
  * @param rep number of repetitions to be performed
  * @return root mean square error of the compared node coordinates
*/
double testElementMapping(mesh *globalMesh, MeshMapping ***mapping, index procIndX,
			  index procIndY, int rep){
    double rmseTotal = 0;
    index nLocalElem = mapping[procIndX][procIndY]->localMesh->nelem;
    
    index randVal, glblElemInd;
    double localCoords[6], globalCoords[6], rmseCurrent;

    // perofrm specified number of repetitions
    for (int i=0; i<rep; ++i){
        // randomly select an index of a local vertex
        randVal = rand() % nLocalElem;
	// get local node coords matching the element with local index randVal
	getLocalNodeCoordsForElement(localCoords, mapping, randVal, procIndX,
				    procIndY);
	// get global element index 
	glblElemInd = mapping[procIndX][procIndY]->elemL2G[randVal];
	// get global node coords matching the element with mapped index glblElemInd
	getGlobalNodeCoordsForElement(globalCoords, globalMesh, glblElemInd,
				      procIndX, procIndY);

	// determine intermediary test result and total RMSE
	rmseCurrent = computeRMSE(globalCoords, localCoords, 6);
	rmseTotal = sqrt(pow(rmseTotal,2.0)+pow(rmseCurrent,2.0));

	DEBUG_PRINT("\nDEBUG ELEMENT MAPPING\n","");
	DEBUG_PRINT("ref:[(%lf,%lf)(%lf,%lf)(%lf,%lf)] test:[(%lf,%lf)(%lf,%lf)(%lf,%lf)] rmse:%lf\n",
		    globalCoords[0],globalCoords[1],globalCoords[2],globalCoords[3],
		    globalCoords[4],globalCoords[5],localCoords[0],localCoords[1],
		    localCoords[2],localCoords[3],localCoords[4],localCoords[5],
		    rmseCurrent);
    }

    return rmseTotal;
}



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
double testVertexMapping(mesh *globalMesh, MeshMapping ***mapping, index procIndX,
                       index procIndY, int rep) {
    
    double rmseTotal = 0;
    index nLocalVert = mapping[procIndX][procIndY]->localMesh->ncoord;

    index randVal;
    double localCoords[2], globalCoords[2], rmseCurrent;
    index mapIndGlbl;

    // perform specified number of repetitions
    for (int i = 0; i < rep; ++i) {
        // randomly select an index of a local vertex
        randVal = rand() % nLocalVert;
        // Get coords of randomly selected local vertex
        localCoords[0] =
            mapping[procIndX][procIndY]->localMesh->coord[randVal * 2];
        localCoords[1] =
            mapping[procIndX][procIndY]->localMesh->coord[randVal * 2 + 1];
        // Get global index through mapping
        mapIndGlbl = mapping[procIndX][procIndY]->vertexL2G[randVal];
        // Get global coordinates
        globalCoords[0] = globalMesh->coord[mapIndGlbl * 2];
        globalCoords[1] = globalMesh->coord[mapIndGlbl * 2 + 1];

	// Determine an intermediary test result and the total RMSE
	rmseCurrent = computeRMSE(globalCoords,localCoords,2); 
	rmseTotal = sqrt(pow(rmseTotal,2.0)+pow(rmseCurrent,2.0));

	DEBUG_PRINT("\nDEBUG VERTEX MAPPING\n","");
	DEBUG_PRINT("ref: (%lf,%lf) tst: (%lf,%lf) err: %lf\n",
		    localCoords[0], localCoords[1], globalCoords[0],
		    globalCoords[1], rmseCurrent);
    }

    return rmseTotal;
}


/**
  * @brief Test for the mapping of vertices, elements and boundaries from the
           local meshes to the global one. Tests are performed for the local
           mesh of each process in the MPI grid
  * @param globalMesh the global unrefined mesh
  * @param gridDims dimensions of the MPI grid
*/
bool testMapping(mesh *globalMesh, int *gridDims) {
    double epsilon = pow(10.0,-15.0);
    bool testResult = 1;

    // initial refinement of the mesh to make the split possible
    index initial_Refines = HPC_MAX(gridDims[0], gridDims[1]);
    mesh *mesh_refined = mesh_initRefinement(globalMesh, initial_Refines);
    // spliet the mesh and thus create a mapping
    MeshMapping ***testMapping = mesh_split(mesh_refined, gridDims);

    printf("%-5s|%-9s|%-9s|%-10s|%-4s\n", "(x,y)","rmse vert","rmse elem",
	    "rmse bndry", "pass");
    printf("-----------------------------------------\n");
    
    // Test the mapping of the local mesh of each process in the MPI grid
    double rmseVertex, rmseElement, rmseBndry;
    bool resultCurrent;

    for (int i = 0; i < gridDims[0]; ++i) {
        for (int j = 0; j < gridDims[1]; ++j) {
	    // Perform the test for the vertices and determine RMSE and itermediary
	    // test result
	    rmseVertex = testVertexMapping(mesh_refined, testMapping, i, j,
			  N_TESTRUNS_MAPPING);
	    resultCurrent = (rmseVertex <= epsilon)?1:0;

	    // Perform the test for the elements
	    rmseElement = testElementMapping(mesh_refined, testMapping, i,j, N_TESTRUNS_MAPPING);
	    resultCurrent = (resultCurrent ==1 && rmseElement <= epsilon)?1:0;

	    // Perform the test for the boundary edges
	    rmseBndry = testBoundaryMapping(mesh_refined, testMapping, i, j, N_TESTRUNS_MAPPING);
	    resultCurrent = (resultCurrent == 1 && rmseBndry <= epsilon)?1:0;

	    // print results
	    printf("(%d,%d)|%-9lf|%-9lf|%-10lf|%-4s\n",i,j,rmseVertex,rmseElement,
		   rmseBndry, resultCurrent?"true":"FALSE");
            printf("-----------------------------------------\n");

	    // update total test result
	    if (testResult && !resultCurrent){
		testResult = resultCurrent;
	    }
        }
    }

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
    int dims[2] = {GRID_N, GRID_M};

    resultCurrent = testNodePartitioning(m1, dims, 8);

    // rewrite total test result only if the test didn't already fail
    if (resultTotal) resultTotal = resultCurrent;

    printf("\n--- Test 2: Mapping ---\n");
    resultCurrent = testMapping(m1, dims);

//    if (resultTotal) resultTotal = resultCurrent;

    printf("\n--- Total Result: %-4s---\n", resultCurrent ? "PASS" : "FAIL");

    printf("=== End test_mesh_split ===\n");
    mesh_free(m1);
}
