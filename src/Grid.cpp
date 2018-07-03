#include <vector>
#include <iostream>

#include "general.h"
#include "Grid.h"

using namespace std;

/*
 * Place a grid covering all the simulation box. Each task keeps track of haloes within some sub-volumes, each grid node 
 * keeps track of all the haloes within the node volume. Haloes are assigned to their nearest grid point.
 * Thus a task needs to activate a sendrecv request to the neighbouring subvolumes/tasks to exchange halos in the buffer volume.
 * The grid is more easily stored as a vector to easily perform an MPI_Gather.
 */

Grid::Grid(){};


Grid::~Grid()
{
	if (locTask == 0)
		cout << "Clearing grid..." << endl;

	taskOnGridNode.clear();
	taskOnGridNode.shrink_to_fit();

	for (int iN = 0; iN < nNodes; iN++)
	{
		haloOnGridNode[iN].clear();
		haloOnGridNode[iN].shrink_to_fit();
	}

	haloOnGridNode.clear();
	haloOnGridNode.shrink_to_fit();

	if (locTask == 0)
		cout << "Done." << endl;
};


void Grid::Init(int n, float size)
{
	N = n; boxSize = size;
	cellSize = (boxSize/n);	// The size of the cell is divided by the total N

	if (locTask == 0)
		cout << "Initialized grid with N=" << N << " nodes in a box of " << boxSize << " kpc/h" << endl; 

	nNodes = N * N * N;

	/* This N^3 vector matrix keeps track on which task holds which part of the grid */
	taskOnGridNode.resize(nNodes);

	/* This vector holds a list of all halos within a given grid cell */
	haloOnGridNode.resize(nNodes);
};


int Grid::Index(int i, int j, int k)
{
	int ii, jj, kk;

	// Add periodic boundary conditions:
	ii = i % N; 	jj = j % N; 	kk = k % N;

	return (ii + N * jj + N * N * kk);
};


int * Grid::GridCoord(float *X)
{
	int * iX; iX = new int[3]; 
	
	for (int i = 0; i<3; i++)
		iX[i] = floor(X[i] / cellSize);

	return iX;
};


/* This function specifies on which task is located each node of the global grid */	// DEPRECATED
void Grid::FindPatchOnTask()
{

/*
	int *ixMax, *ixMin;
	ixMax = new int[3];	ixMin = new int[3];	
	
	ixMax = GridCoord(locXmax);
	ixMin = GridCoord(locXmin);

	for (int i = ixMin[0]; i < ixMax[0]; i++)
		for (int j = ixMin[1]; j < ixMax[1]; j++)
			for (int k = ixMin[2]; k < ixMax[2]; k++)
				taskOnGridNode[Index(i, j, k)] = locTask + 1;	// distinguish from the "null" nodes

	free(ixMax); free(ixMin);

*/
};


void Grid::AssignToGrid(float *X, int index)
{
	int *iX, thisNode = 0;
	iX = new int[3];	

	iX = GridCoord(X);
	thisNode = Index(iX[0], iX[1], iX[2]);
	haloOnGridNode[thisNode].push_back(index);
	taskOnGridNode[thisNode] = locTask + 1;	// Add one to distinguish from empty node!

//	if (index < 50) 	// Sanity check
//		printf("%d) Halo=%d grid=(%d, %d, %d) x=(%.2f, %.2f, %.2f) node=%d\n", 
//			locTask, index, iX[0], iX[1], iX[2], X[0], X[1], X[2], thisNode);
	free(iX);
};


void Grid::RecvBufferNodes()
{

#ifdef TEST_BLOCK
	// On each task, we identify the edge of the subbox and identify the additional nodes needed 

	// Create an array with all the tasks, which contains the list of nodes to be requested from each task
	// recvNodes[0] = [18922929, 10238933, 192929] ---> these are the nodes that locTask will request from task 0
	// do a push_back() for each new node encountered

	// PSEUDO-CODE FIXME TODO
	vector <vector <int>> recvNodes;
	recvNodes.resize(totTask);

	nCellsBuffer = VmaxTot * timeStep / cellSize;

	xMinBuff = xMin - nCellsBuff;
	for (ix = xMinBuff; ix < xMin; ix ++)
	// Take a slab - fix the ix and then loop over all 
	//r (iy...)
	thisNode = Index(ix, iy, iz);
	thisTask = taskOnGridNode[thisNode];	
	recvNodes[thisTask].push_back(thisNode);
	
//		for (iy = xMinBuff; ix < xMin; ix ++)
 
#endif
};


void Grid::Info()
{	
	int *ixMax, *ixMin;
	ixMax = new int[3];	ixMin = new int[3];	
/*
	for (int i = ixMin[0]; i < ixMax[0]; i++)
		for (int j = ixMin[1]; j < ixMax[1]; j++)
			for (int k = ixMin[2]; k < ixMax[2]; k++)
				cout << Index(i, j, k) << " " << taskOnGridNode[Index(i, j, k)] << endl; 
*/
};


