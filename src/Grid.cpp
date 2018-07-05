#include <vector>
#include <iostream>
#include <algorithm>

#include "general.h"
#include "Grid.h"

using namespace std;

/*
 * Place a grid covering all the simulation box. Each task keeps track of haloes within some sub-volumes, each grid node 
 * keeps track of all the haloes within the node volume. Haloes are assigned to their nearest grid point.
 * Thus a task needs to activate a sendrecv request to the neighbouring subvolumes/tasks to exchange halos in the buffer volume.
 * The grid is more easily stored as a vector to easily perform an MPI_Gather.
 */

Grid::Grid()
{
};


Grid::~Grid()
{
	//if (locTask == 0)
	//	cout << "Clearing grid..." << endl;

	taskOnGridNode.clear();
	taskOnGridNode.shrink_to_fit();

	for (int iN = 0; iN < nNodes; iN++)
	{
		haloOnGridNode[iN].clear();
		haloOnGridNode[iN].shrink_to_fit();
	}

	haloOnGridNode.clear();
	haloOnGridNode.shrink_to_fit();

	//if (locTask == 0)
	//	cout << "Done." << endl;
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


int * Grid::Index2Grid(int index)
{
	int *iX; iX = new int[3];

	iX[2] = floor(index / (N * N));
	iX[1] = floor((index - iX[2] * N * N) / N);
	iX[0] = floor((index - iX[2] * N * N - iX[1] * N));

	return iX;	
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


void Grid::FindNearbyNodes(int index, int nCells)
{
	int thisIndex, thisTask;
	int *iX, *jX; jX = new int[3]; 

	iX = Index2Grid(index);

	/* Periodic boundary conditions are taken care of by the Index() function, no need to implement them here */
	for (int i = -nCells; i < nCells+1; i++)
	{
		//jX[0] = iX[0] + i % N;
		jX[0] = iX[0] + i;

		for (int j = -nCells; j < nCells+1; j++)
		{
			jX[1] = iX[1] + j;

			for (int k = -nCells; k < nCells+1; k++)
			{ 
				jX[2] = iX[2] + k;
				
				thisIndex = Index(jX[0], jX[1], jX[2]);

				/* There might be nodes shared among several tasks, so we have to loop here */
				for (int l = 0; l < globalTaskOnGridNode[thisIndex].size(); l++)
				{
					// Tasks on the grid are named with a +1 to avoid confusion with null nodes
					thisTask = globalTaskOnGridNode[thisIndex][l] -1;
					//cout << thisTask << " " << thisIndex << endl;

					// Only add the task to the communication buffer if the node is 
					if (thisTask != locTask && thisTask > -1)
					{
						buffNodes[thisTask].push_back(thisIndex);
						//cout << locTask << " needs " << thisIndex << " from task " 
						//	<< thisTask << " " << endl;
					}
				}
			}
		}
	}
};


void Grid::SortLocNodes()
{
	sort(locNodes.begin(), locNodes.end() );
	locNodes.erase(unique(locNodes.begin(), locNodes.end()), locNodes.end());
/*
	for (int i = 0; i < locNodes.size(); i++)
	{
		int *iX = Index2Grid(locNodes[i]);
		cout << "Task " << locTask << " has node: " << locNodes[i] <<  " GridX: " << iX[0] << " " << iX[1] << " " << iX[2] << endl;
	}
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

	locNodes.push_back(thisNode);

//	if (index < 50) 	// Sanity check
//		printf("%d) Halo=%d grid=(%d, %d, %d) x=(%.2f, %.2f, %.2f) node=%d\n", 
//			locTask, index, iX[0], iX[1], iX[2], X[0], X[1], X[2], thisNode);
	free(iX);
};


/* 
 * This function identifies the nodes of the buffer region for each task.
 * For each halo it determines the nodes allocated to other tasks and stores them into a list.
 * It loops on all the nodes already allocated and identifies all the neighbouring nodes within 
 * a radius of maxBufferThick size.
 */
void Grid::FindBufferNodes()
{
	int nCells = 1;

	if (buffNodes.size() == 0)
		buffNodes.resize(totTask);
	//cout << locTask << "- size " << buffNodes.size() << endl;;

	// Do a loop on all the nodes contained in this task to find out which nodes need to be communicated
	for (int i = 0; i < locNodes.size(); i++)
		FindNearbyNodes(locNodes[i], nCells);		
	
	/*	Clean Buffer Nodes	*/
	for (int i = 0; i < totTask; i++)
	{
		sort(buffNodes[i].begin(), buffNodes[i].end());
		buffNodes[i].erase(unique(buffNodes[i].begin(), buffNodes[i].end()), buffNodes[i].end());
		//cout << locTask << " needs " << buffNodes[i].size() << " nodes from task=" << i << endl;
	}

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


