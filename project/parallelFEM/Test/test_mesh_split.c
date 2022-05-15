#include <assert.h>

#include "hpc.h"

#define GRID_M 2
#define GRID_N 2

void getTargetNodeCounts(index *buffer,mesh *globalMesh, 
			      index gridDimX,index gridDimY){
    for (index i=0; i<gridDimY; ++i){

	for (index j=0; j<gridDimX; ++j){
	    buffer[i*gridDimX+j] = getSliceSize((index)sqrt((double)globalMesh->ncoord - 1),
					    gridDimY, i)+1;
	    buffer[i*gridDimX+j] *= getSliceSize((index)sqrt((double)globalMesh->ncoord - 1),
					    gridDimX, j)+1;
	}
    }
}

void collectLocalNodeCounts(index *buffer, MeshMapping ***mapping,
			    index gridDimX,index gridDimY){
    for (int i=0; i<gridDimY; ++i){
	for (int j=0; j<gridDimX; ++j){
	    buffer[i*gridDimX+j] = mapping[i][j]->localMesh->ncoord;
	}
    }
}

void getTargetBndryCounts(index *buffer, mesh *mesh, index gridDimX, index gridDimY){
    index nEdgesX, nEdgesY; 

    for (index i=0; i<gridDimY; ++i){
	nEdgesY = getSliceSize((index)sqrt((double)mesh->ncoord - 1),
			       gridDimY, i);

	for (index j=0; j<gridDimX; ++j){
	    nEdgesX = getSliceSize((index)sqrt((double)mesh->ncoord - 1),
				    gridDimX, j);

	    buffer[i*gridDimX+j] = 0;

	    if (i==0 || i==gridDimY-1){
		buffer[i*gridDimX+j] = nEdgesX;
	    } 
	    if (j==0 || j==gridDimX-1){
		buffer[i*gridDimX+j] += nEdgesY;
	   }
	}

    }
}

void collectLocalBndryCounts(index *buffer, MeshMapping ***mapping,
			     index gridDimX, index gridDimY){
     for (int i=0; i<gridDimY; ++i){
	for (int j=0; j<gridDimX; ++j){
	    buffer[i*gridDimX+j] = mapping[i][j]->localMesh->nbdry;
	}
    }
}


mesh *initRefinement(mesh *globalMesh, index *gridDims){
    index n = HPC_MAX(gridDims[0],gridDims[1]);
    mesh *mesh_current;
    mesh *mesh_previous;

    mesh_current = mesh_refine(globalMesh);
    mesh_getEdge2no(mesh_current->nelem, mesh_current->elem, &mesh_current->nedges,
		    &mesh_current->edge2no);
    mesh_current->fixed = mesh_getFixed(mesh_current->ncoord, mesh_current->bdry,
					mesh_current->nbdry, &mesh_current->nfixed);

    for (int i=1; i<n; ++i){
	mesh_previous = mesh_current;
	mesh_current = mesh_refine(mesh_current);
	mesh_getEdge2no(mesh_current->nelem, mesh_current->elem, &mesh_current->nedges,
		        &mesh_current->edge2no);
	mesh_current->fixed = mesh_getFixed(mesh_current->ncoord, mesh_current->bdry,
					    mesh_current->nbdry, &mesh_current->nfixed);
	mesh_free(mesh_previous);
    }

    return mesh_current;
}



void testNodePartitioning(mesh *globalMesh, index *gridDims, int numRefines){
    mesh *mesh_current;
    mesh *mesh_previous;

    // initial refinement of the mesh to make the mapping possible according to
    // the mpi grid
    mesh_current = initRefinement(globalMesh, gridDims);
    index initialRefines = HPC_MAX(gridDims[0],gridDims[1]);

    printf("\n%-4s|%-5s|%-14s|%-14s|%-8s|%-8s|%-5s\n",
	    "ref","(i,j)","cnt nodes ref","cnt nodes tst","bdry ref","bdry tst","pass");
    printf("---------------------------------------------------------------\n");

    // Refine mesh iteratively and collect target and local node counts
    index refNodeCounts[gridDims[0]*gridDims[1]];
    index testNodeCounts[gridDims[0]*gridDims[1]];
    index refBndryCounts[gridDims[0]*gridDims[1]];
    index tstBndryCounts[gridDims[0]*gridDims[1]];
    MeshMapping ***testMapping;

    for (int i=initialRefines; i<numRefines; ++i){
	// create mapping
	testMapping = mesh_split(mesh_current, gridDims);

	// Get raference node counts and node counts of local meshes
	getTargetNodeCounts(refNodeCounts, mesh_current, gridDims[0], gridDims[1]);
	collectLocalNodeCounts(testNodeCounts, testMapping, gridDims[0],gridDims[1]);

	// Get refenerence and local bndry edge counts 
	getTargetBndryCounts(refBndryCounts, mesh_current, gridDims[0], gridDims[1]);
	collectLocalBndryCounts(tstBndryCounts, testMapping, gridDims[0], gridDims[1]); 

	// Compare and print result
	index refCount, tstCount, refBndryCnts, tstBndryCnts;
	for (int k=0; k<gridDims[1]; ++k){
	    for (int l=0; l<gridDims[0]; ++l){
		refCount = refNodeCounts[k*gridDims[0]+l];
		tstCount = testNodeCounts[k*gridDims[0]+l];
		refBndryCnts = refBndryCounts[k*gridDims[0]+l];
		tstBndryCnts = tstBndryCounts[k*gridDims[0]+l];
		printf("%-4td|(%td,%td)|%-14td|%-14td|%-8td|%-8td|%-5s\n",
		      i, k, l, refCount, tstCount, refBndryCnts, tstBndryCnts,
		      (refCount==tstCount && refBndryCnts == tstBndryCnts)?"true":"false");
	    }
	}
	printf("---------------------------------------------------------------\n");
	// next refinement
	mesh_previous = mesh_current;
	mesh_current = mesh_refine(mesh_previous);
	mesh_getEdge2no(mesh_current->nelem, mesh_current->elem, &mesh_current->nedges,
		    &mesh_current->edge2no);
	mesh_current->fixed = mesh_getFixed(mesh_current->ncoord, mesh_current->bdry,
		    mesh_current->nbdry, &mesh_current->nfixed);
	
	// delete old mapping
	delete2DMeshMapping(testMapping, gridDims[0]);

	// free old mesh
	mesh_free(mesh_previous);

    }

    // free final mesh
    mesh_free(mesh_current);

}

int main() {
    printf("\n=== Start test_mesh_split ===\n");

    char fname[32] = "../Problem/problem1";
    mesh *m1 = mesh_load(fname);
    mesh_getEdge2no(m1->nelem, m1->elem, &m1->nedges, &m1->edge2no);
    m1->fixed = mesh_getFixed(m1->ncoord, m1->bdry, m1->nbdry, &m1->nfixed);
   
    printf("\n--- Test 1: Partitioning of Nodes ---");
    index dims[2] = {2,4};

    testNodePartitioning(m1, dims, 8);



    printf("Okay\n");
    printf("=== End test_mesh_split ===\n");

    mesh_free(m1);
}
