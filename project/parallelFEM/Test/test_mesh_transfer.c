#include "hpc.h"
#include <math.h>
#include <string.h>

#ifndef MPI_GRID_X 
#define MPI_GRID_X  2
#endif

#ifndef MPI_GRID_Y
#define MPI_GRID_Y  2
#endif

#ifndef NUM_MESH_REFINES
#define NUM_MESH_REFINES 2
#endif

/**
 * @brief function to write the mesh/mapping metadata for the mapping of the local
 mesh of a process from the global mapping to a buffer
 * @param buffer pointer to array of type index with length 6
 * @param globalMapping pointer to two dimensional array of MeshMapping structs
 i.e. the global mapping conatining the mapping for all processes
 * @param dim coordinates of the process in the MPI grid
 */
void getReferenceMetaData(index *buffer, MeshMapping ***globalMapping, int *dims){
    buffer[0] = globalMapping[dims[0]][dims[1]]->localMesh->ncoord;
    buffer[1] = globalMapping[dims[0]][dims[1]]->localMesh->nelem;
    buffer[2] = globalMapping[dims[0]][dims[1]]->localMesh->nedges;
    buffer[3] = globalMapping[dims[0]][dims[1]]->localMesh->nbdry;
    buffer[4] = globalMapping[dims[0]][dims[1]]->localMesh->nfixed;
    buffer[5] = globalMapping[0][0]->globalNcoord; 
}

/** 
 * @brief function to compute the root mean square error of test data with type 
 index
 * @param refData array of reference data
 * @param tstData array of test data
 * @param numElements number of entries in the arrays
 * @return root mean square error of the test data
 */
double index_computeRMSE(index *refData, index *tstData, int numElements){
    double rmse=0;
    for (int i=0; i<numElements; ++i){
	rmse += pow((double)(refData[i]-tstData[i]),2.0);
    }
    rmse = sqrt(rmse);
    return rmse;
}

/**
 * @brief function for testing the transfer of the mesh metadata. For the test,
 the test and reference data are first collected from the global and 
 local data on the respective process. The testdata is are then communicated
 to the root process and the rmse of the data is evaluated.
 * @param globalMesh the global mesh struct
 * @param globalMapping pointer to the global mapping, containig the MeshMapping structs
 for all processes. Can be Nullpointer for nonroot processes
 * @param localMapping pointer to the local mapping struct of the respective nonroot process
 * @param grid the MPI grid
 * @return boolean value stating the sucess of the test
 */
bool testMeshMetadata(mesh *globalMesh, MeshMapping ***globalMapping,
	MeshMapping *localMapping, MPI_Comm grid){
    int rank; MPI_Comm_rank(grid, &rank);
    int nof_processes; MPI_Comm_size(grid, &nof_processes);
    bool testResult = 1;
    char testRestultString[12] = "FALSE";
    double epsilon = pow(10,-15);

    if (rank==0){
	printf("%-4s|(x,y)|%-8s|%-8s|%-8s|%-8s|%-8s|%-11s|%-5s|%-5s\n","","ncoord","nelem","nedges",
		"nbdry","nfixed","glbl ncoord","RMSE","pass");
	printf("----------------------------------------------------------------------------------\n");

	index refMetadata[6], tstMetadata[6];
	int procCoords[2];
	double rmse;
	bool resultCurrent;

	// iterate over all nonroot processes
	for (int r=0; r<nof_processes; ++r){
	    // get coords of processes with rank r in MPI grid
	    MPI_Cart_coords(grid, r, 2, procCoords);

	    // get reference metadata of the process
	    getReferenceMetaData(refMetadata, globalMapping, procCoords);
	    printf("%-4s|(%d,%d)|%-8td|%-8td|%-8td|%-8td|%-8td|%-11td|%-5d|%-5s\n",
		    "ref",procCoords[0],procCoords[1],refMetadata[0],refMetadata[1],
		    refMetadata[2], refMetadata[3],refMetadata[4],refMetadata[5],0," ");

	    // get metadata from the mapping of current nonroot process by communication
	    // and compute RMSE
	    MPI_Status status;
	    if (rank==0){
		tstMetadata[0] = localMapping->localMesh->ncoord;
		tstMetadata[1] = localMapping->localMesh->nelem;
		tstMetadata[2] = localMapping->localMesh->nedges;
		tstMetadata[3] = localMapping->localMesh->nbdry;
		tstMetadata[4] = localMapping->localMesh->nfixed;
		tstMetadata[5] = localMapping->globalNcoord;
	    } else{
		MPI_Recv(tstMetadata, 6, MPI_LONG_LONG, r, 1, grid, &status);
	    }
	    rmse = index_computeRMSE(refMetadata, tstMetadata, 6);
	    if (rmse<epsilon){
		strcpy(testRestultString,"true");
		testResult = 0;
	    } else {
		strcpy(testRestultString,"FALSE");
		testResult	= 1;
	    }
	    printf("%-4s|(%d,%d)|%-8td|%-8td|%-8td|%-8td|%-8td|%-11td|%-1.3lf|%-5s\n",
		    "tst",procCoords[0],procCoords[1],tstMetadata[0],tstMetadata[1],
		    tstMetadata[2], tstMetadata[3],tstMetadata[4],tstMetadata[5],
		    rmse,testRestultString);
	}
    }else{
	// Collect metadata of the local mapping
	index sendData[6] = {localMapping->localMesh->ncoord,
	    localMapping->localMesh->nelem,
	    localMapping->localMesh->nedges,
	    localMapping->localMesh->nbdry,
	    localMapping->localMesh->nfixed,
	    localMapping->globalNcoord};
	// Send the metadata of the local mapping
	MPI_Send(sendData, 6, MPI_LONG_LONG, 0, 1, grid);
    }
    return testResult;
}

/**
 * @brief function to compute the RMSE of testdata of type double with respect
 to reference data 
 * @param refData array of reference data
 * @param tstData array of test data
 * @param numElements number of elements in the data arrays
 * @return RMSE of the test data
 */
double double_computeRMSE(double *refData, double *tstData, int numElements){
    double rmse=0;
    for (int i=0; i<numElements; ++i){
	rmse += pow(refData[i]-tstData[i],2.0);
    }
    rmse = sqrt(rmse);
    return rmse;
}

double print_coords(double *coords, index ncoord){
    for (index i=0; i<ncoord; ++i){
	printf("(%lf,%lf)\n",coords[i*2],coords[i*2+1]);
    }
}


/**
 * @brief function for testing the transfer of data for the local meshes from
 root to the nonroot processes. For the tests, each data array of the
 local mesh, stored on the respective nonroot process is compared to
 the original data on the root process.
 * @param globalMapping pointer to the global mapping, that contains the mesh/mapping
 data of all processes. Can be Nullpointer for nonroot processes
 * @param localMapping pointer to the local mapping struct of the respective nonroot process
 * @param globalMesh pointer to the global mesh struct
 * @param grid the MPI grid
 * @return boolean value stating the sucess of the test
 */
bool testMeshData(MeshMapping ***globalMapping, MeshMapping *localMapping,
	mesh *globalMesh, MPI_Comm grid){
    int rank; MPI_Comm_rank(grid, &rank);
    int nof_processes; MPI_Comm_size(grid, &nof_processes);
    bool testResult = 1;
    double epsilon = pow(10,-15);

    if (rank==0){
	// allocate sufficient memory for the local Data
	index *testElemBdryEdgeData = malloc(7*globalMesh->nelem * sizeof(index));
	double *testCoordData = malloc(2*globalMesh->ncoord * sizeof(index));
	if (!testElemBdryEdgeData && !testCoordData){
	    printf("\nMemory allocation for test data failed! ABORT!\n");
	    abort();
	}

	printf("(x,y)|%-10s|%10s|%10s|%12s|%10s|%-5s","RMSE coord","RMSE elem",
		"RMSE bdry","RMSE edge2no","RMSE fixed","pass\n");
	printf("-------------------------------------------------------------------\n");

	int pCoords[2];
	double rmseCoord, rmseElem, rmseBdry, rmseEdge2no, rmseFixed;
	bool resultCurrent;
	index nofData;
	MPI_Status status;

	// Perform tests for all processes iteratively 
	for (int r=0; r<nof_processes; ++r){
	    // get coords of processes with rank r in MPI grid
	    MPI_Cart_coords(grid, r, 2, pCoords);  


	    // Test for coordinate data

	    // Set the value in the memory to zero, such that smaller sizes of
	    // the data do not matter for the rmse
	    nofData = 2*globalMapping[pCoords[0]][pCoords[1]]->localMesh->ncoord;
	    memset(testCoordData, 0.0, nofData*sizeof(double));
	    /* get the test data by communicating or copying in case of root,
	       and compute RMSE*/
	    if (r ==0){
		memcpy(testCoordData, localMapping->localMesh->coord,
			nofData*sizeof(double));
	    } else{
		MPI_Recv(testCoordData,nofData, MPI_DOUBLE, r, 2, grid, &status);
	    }
	    rmseCoord = double_computeRMSE(globalMapping[pCoords[0]][pCoords[1]]
		    ->localMesh->coord, testCoordData, nofData);

	    // Test for element data
	    nofData = 7*globalMapping[pCoords[0]][pCoords[1]]->localMesh->nelem;
	    memset(testElemBdryEdgeData, 0, nofData*sizeof(index));
	    if (r==0){
		memcpy(testElemBdryEdgeData, localMapping->localMesh->elem,
			nofData*sizeof(index));
	    } else{
		MPI_Recv(testElemBdryEdgeData, nofData, MPI_LONG_LONG, r, 3, grid,
		    &status);
	    }
	    rmseElem = index_computeRMSE(globalMapping[pCoords[0]][pCoords[1]]
		    ->localMesh->elem, testElemBdryEdgeData,
		    nofData);

	    // Test for boundary data
	    nofData = 4*globalMapping[pCoords[0]][pCoords[1]]->localMesh->nbdry;
	    memset(testElemBdryEdgeData, 0, nofData*sizeof(index));
	    if (r==0){
		memcpy(testElemBdryEdgeData, localMapping->localMesh->bdry,
			nofData*sizeof(index));
	    } else{
		MPI_Recv(testElemBdryEdgeData, nofData, MPI_LONG_LONG, r, 4, grid,
			 &status);
	    }
	    rmseBdry = index_computeRMSE(globalMapping[pCoords[0]][pCoords[1]]
		    ->localMesh->bdry, testElemBdryEdgeData,
		    nofData);

	    // Test for edge2no data only if pointer is defined
	    if (globalMapping[pCoords[0]][pCoords[1]]->localMesh->edge2no){
		nofData = 2*globalMapping[pCoords[0]][pCoords[1]]->localMesh->nedges;
		memset(testElemBdryEdgeData, 0, nofData*sizeof(index));
		if (r==0){
		    memcpy(testElemBdryEdgeData, localMapping->localMesh->edge2no,
			    nofData*sizeof(index));
		} else{
		    MPI_Recv(testElemBdryEdgeData, nofData, MPI_LONG_LONG, r, 5, grid,
			    &status);
		}
		rmseEdge2no = index_computeRMSE(globalMapping[pCoords[0]][pCoords[1]]
			->localMesh->edge2no, testElemBdryEdgeData,
			nofData);
	    } else{
		rmseEdge2no = -1;
	    }

	    // test for fixed boudaries only if pointer is defined
	    if (globalMapping[pCoords[0]][pCoords[1]]->localMesh->fixed){
		nofData = 4*globalMapping[pCoords[0]][pCoords[1]]->localMesh->nfixed;
		memset(testElemBdryEdgeData, 0, nofData*sizeof(index));
		if (r==0){
		    memcpy(testElemBdryEdgeData, localMapping->localMesh->fixed,
			    nofData*sizeof(index));
		} else{
		    MPI_Recv(testElemBdryEdgeData, nofData, MPI_LONG_LONG, r, 6, grid,
				&status);
		}
		rmseFixed = index_computeRMSE(globalMapping[pCoords[0]][pCoords[1]]
			->localMesh->fixed, testElemBdryEdgeData,
			nofData);
	    } else{
		rmseFixed = -1;
	    }

	    // determine intermediary result and print RMSEs
	    resultCurrent = (rmseCoord <= epsilon && rmseElem <= epsilon &&
		    rmseBdry <= epsilon && rmseEdge2no <= epsilon &&
		    rmseFixed <= epsilon)?1:0
		; 
	    printf("(%d,%d)|%10lf|%10lf|%10lf|%12lf|%10lf|%-5s\n",pCoords[0],pCoords[1],
		    rmseCoord,rmseElem,rmseBdry,rmseEdge2no,rmseFixed,
		    resultCurrent?"true":"FALSE");

	    // update test result
	    if (testResult){
		testResult = resultCurrent;
	    }

	}

	// free the memory for the test data
	free(testElemBdryEdgeData);
	free(testCoordData);
    } else{
	// Send coordinate data to root
	MPI_Send(localMapping->localMesh->coord,2*localMapping->localMesh->ncoord,
		MPI_DOUBLE, 0, 2, grid);

	// Send element data to root
	MPI_Send(localMapping->localMesh->elem, 7*localMapping->localMesh->nelem,
		MPI_LONG_LONG, 0, 3, grid);

	// Send boundary data to root
	MPI_Send(localMapping->localMesh->bdry,4*localMapping->localMesh->nbdry,
		MPI_LONG_LONG, 0, 4, grid);

	// If existing: Send edge2no data to root
	if (localMapping->localMesh->edge2no){
	    MPI_Send(localMapping->localMesh->edge2no,2* localMapping->localMesh->nedges,
		    MPI_LONG_LONG, 0, 5, grid);
	}

	// If existing: Send fixed boundary data to root
	if (localMapping->localMesh->fixed){
	    MPI_Send(localMapping->localMesh->fixed,4* localMapping->localMesh->nfixed,
		    MPI_LONG_LONG, 0, 6, grid);
	}
    }

    return testResult;
}


/**
 * @brief function to test the transfer of the mapping data for the local mapping
 of a process from its origin on root to the respective processs. For
 the test, the local data are tansferred back from the respective nonroot
 process and compared to the original data on the root process, by computing
 the RMSE.
 * @param globalMapping pointer to the global mapping, containing the mapping data
 for all processes. Can be nullpointer for nonroot processes calling the
 function
 * @param localMapping pointer to local mapping struct of the respective process
 * @param grid the MPI grid
 * @return boolean value stating the sucess of the test
 */
bool testMappingData(MeshMapping ***globalMapping, MeshMapping *localMapping,
	MPI_Comm grid){
    int rank; MPI_Comm_rank(grid, &rank);
    int nof_processes; MPI_Comm_size(grid, &nof_processes);
    bool testResult = 1;
    double epsilon = pow(10,-15);

    if (rank==0){
	// Allocate sufficient memory space for the test data
	index *testMappingData = malloc(localMapping->globalNcoord*sizeof(index));
	if (!testMappingData){
	    printf("\nAllocation of memory for test data failed! ABORT!\n");
	    abort();
	}

	// Header for Test output
	printf("(x,y)|%10s|%10s|%10s|%5s\n","RMSE Vert","RMSE Elem","RMSE Bdry",
		"pass");
	printf("--------------------------------------------\n");


	// Test all processes iteratively
	int pCoords[2];
	double rmseVertex, rmseElem, rmseBdry;
	bool resultCurrent;
	index nofData;
	MPI_Status status;

	for (int r=0; r<nof_processes; ++r){
	    // get coordinates of current process
	    MPI_Cart_coords(grid, r, 2, pCoords);


	    // Test for vertex mapping data
	    nofData = globalMapping[pCoords[0]][pCoords[1]]->localMesh->ncoord;
	    memset(testMappingData, 0, localMapping->globalNcoord*sizeof(index));
	    if (r==0){
		memcpy(testMappingData, localMapping->vertexL2G,
			nofData*sizeof(index));
	    } else{
		MPI_Recv(testMappingData, nofData, MPI_LONG_LONG, r, 7, grid, &status);
	    }
	    rmseVertex = index_computeRMSE(globalMapping[pCoords[0]][pCoords[1]]
		    ->vertexL2G, testMappingData, nofData);

	    // Test for element mapping data
	    nofData = globalMapping[pCoords[0]][pCoords[1]]->localMesh->nelem;
	    memset(testMappingData, 0, localMapping->globalNcoord*sizeof(index));
	    if (r==0){
		memcpy(testMappingData, localMapping->elemL2G,
			nofData*sizeof(index));
	    } else{
		MPI_Recv(testMappingData, nofData, MPI_LONG_LONG, r, 8, grid, &status);
	    }
	    rmseElem = index_computeRMSE(globalMapping[pCoords[0]][pCoords[1]]
		    ->elemL2G, testMappingData, nofData);

	    // Test for boundary mapping data
	    nofData = globalMapping[pCoords[0]][pCoords[1]]->localMesh->nbdry;
	    memset(testMappingData, 0, localMapping->globalNcoord*sizeof(index));
	    if (r==0){
		memcpy(testMappingData, localMapping->bdryL2G,
			nofData*sizeof(index));
	    } else{
		MPI_Recv(testMappingData, nofData, MPI_LONG_LONG, r, 9, grid, &status);
	    }
	    rmseBdry = index_computeRMSE(globalMapping[pCoords[0]][pCoords[1]]
		    ->bdryL2G, testMappingData, nofData);

	    // Determine intermediary result and print RMSEs
	    resultCurrent = (rmseVertex <= epsilon && rmseElem <= epsilon &&
		    rmseBdry <= epsilon)?1:0;

	    printf("(%d,%d)|%10lf|%10lf|%10lf|%5s\n", pCoords[0],pCoords[1],
		    rmseVertex, rmseElem, rmseBdry, resultCurrent?"true":"FALSE");

	    // update total test result
	    if (testResult){
		testResult = resultCurrent;
	    }

	}

    } else{
	// Send mapping data of vertices to root
	MPI_Send(localMapping->vertexL2G, localMapping->localMesh->ncoord,
		MPI_LONG_LONG, 0, 7, grid);

	// Send mapping data of elements to root
	MPI_Send(localMapping->elemL2G, localMapping->localMesh->nelem,
		MPI_LONG_LONG, 0, 8, grid);

	// Send mapping data of boundary edges to root
	MPI_Send(localMapping->bdryL2G, localMapping->localMesh->nbdry,
		MPI_LONG_LONG, 0, 9, grid);
    }

    return testResult;

}


int main(int argc, char**argv){
    MPI_Init(&argc, &argv);

    // Initialize MPI gird Communication
    int nof_processes; MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int dims[2] = {MPI_GRID_X,MPI_GRID_Y};
    int periods[2] = {false,false};

    bool resultTotal = 1;

    char fname[32] = "../Problem/problem1";

    MPI_Dims_create(nof_processes, 2, dims);
    MPI_Comm grid;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &grid);
    MPI_Comm_rank(grid, &rank);


    MeshMapping ***mapping;
    mesh *m1;
    mesh *meshRefined;
    MeshMapping *localMapping;

    // root: load and refine mesh, send local meshes
    if (rank==0){
	printf("\n=== Start test_mesh_transfer ==\n");	

	// Load basic mesh
	mesh *m1 = mesh_load(fname);
	mesh_getEdge2no(m1->nelem, m1->elem, &m1->nedges, &m1->edge2no);
	m1->fixed = mesh_getFixed(m1->ncoord, m1->bdry, m1->nbdry, &m1->nfixed);

	// initial refinement
	meshRefined = mesh_initRefinement(m1, NUM_MESH_REFINES); 

	// create mapping
	mapping = mesh_split(meshRefined, dims);
    }

    // transfer local meshes
    localMapping = mesh_transfer(mapping, grid);

    // Root: do the testing
    if (rank==0){
	bool resultCurrent;
	bool result = 1;
	index meshDim=(index)sqrt((double)meshRefined->ncoord);
	printf("\nNumber of mesh refinements: %d\n",NUM_MESH_REFINES);
	printf("Mesh dimesnion: (%td,%td)\n",meshDim,meshDim);


	// start Test for Metadata
	printf("\n--- Test 1: Metadata ---\n");	
	resultCurrent = testMeshMetadata(meshRefined, mapping, localMapping,  grid);
	// update total result only if a test hasn't already failed
	result = result?resultCurrent:result;

	// start Test for Coorinate data
	printf("\n--- Test 2: Meshdata ---\n");
	resultCurrent = testMeshData(mapping, localMapping, meshRefined, grid);
	result = result?resultCurrent:result;

	// start Test for Mapping data
	printf("\n--- Test 3: Mapping data ---\n");
	resultCurrent = testMappingData(mapping,localMapping,grid);
	result = result?resultCurrent:result;

	printf("\n--- Test result: %s\n",result==0?"PASS":"FAILED");

	// Nonroot: call the test functions to send data for testing
    } else{
	testMeshMetadata(meshRefined, mapping, localMapping,  grid);

	testMeshData(mapping, localMapping, meshRefined, grid);

	testMappingData(mapping,localMapping,grid);
    }


    // If rank=0 free the mapping, the refined mesh and the base mesh
    if (rank ==0){
	delete2DMeshMapping(mapping, dims[0]);
	mesh_free(meshRefined);
	// mesh_free(m1);

	// if rank!=0 free the local mesh
    } else{
	if (localMapping != NULL){
	    mesh_free(localMapping->localMesh);
	    deleteMeshMapping(localMapping);
	}
    }

    DEBUG_PRINT("\nProcess with rank %d ending Test!\n",rank);    

    if (rank==0){
	printf("\n=== End test_mesh_transfer ===\n");
    }

    MPI_Finalize();
}
