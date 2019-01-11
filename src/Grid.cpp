/*
 *   METROC++: MErger TRees On C++, a scalable code for the computation of merger trees in cosmological simulations.
 *   Copyright (C) Edoardo Carlesi 2018-2019
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/*
 * Grid.cpp:
 * Functionalites for the grid covering the full simulation box.
 * Halos are assigned to their nearest grid point (NGP) and nodes are communicated all across the tasks.
 * In this way, each node knows to which task (and to which task) it needs to send (and receive) halos in the 
 * buffer region.
 */

#include <vector>
#include <iostream>
#include <algorithm>

#include "utils.h"
#include "global_vars.h"
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
	Clean();
};


void Grid::Clean()
{
	if (locNodes.size() > 0)
	{
		locNodes.clear();
		locNodes.shrink_to_fit();
	}

	if (taskOnGridNode.size() > 0)
	{
		taskOnGridNode.clear();
		taskOnGridNode.shrink_to_fit();
	}

	if (globalTaskOnGridNode.size() > 0)
	{

		for (int iN = 0; iN < globalTaskOnGridNode.size(); iN++)
		{
			globalTaskOnGridNode[iN].clear();
			globalTaskOnGridNode[iN].shrink_to_fit();
		}

		globalTaskOnGridNode.clear();
		globalTaskOnGridNode.shrink_to_fit();
	}

	if (haloOnGridNode.size() > 0)
	{
		for (int iN = 0; iN < haloOnGridNode.size(); iN++)
		{
			haloOnGridNode[iN].clear();
			haloOnGridNode[iN].shrink_to_fit();
		}

		haloOnGridNode.clear();
		haloOnGridNode.shrink_to_fit();
	}

	if (buffOnGridNode.size() > 0)
	{
		for (int iN = 0; iN < buffOnGridNode.size(); iN++)
		{
			buffOnGridNode[iN].clear();
			buffOnGridNode[iN].shrink_to_fit();
		}

		buffOnGridNode.clear();
		buffOnGridNode.shrink_to_fit();
	}

	if (buffNodes.size() > 0)
	{
		for (int iN = 0; iN < buffNodes.size(); iN++)
		{
			buffNodes[iN].clear();
			buffNodes[iN].shrink_to_fit();
		}

		buffNodes.clear();
		buffNodes.shrink_to_fit();	
	}
};


void Grid::Init(int n, float size)
{
	N = n; boxSize = size;
	cellSize = (boxSize/n);	// The size of the cell is divided by the total N

#ifdef VERBOSE
	if (locTask == 0)
		cout << "Initialized grid with N=" << N << " nodes in a box of " << boxSize << " kpc/h" << endl; 
#endif

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
	ii = (i + N) % N; 	jj = (j + N) % N; 	kk = (k + N) % N;

	return (ii + N * jj + N * N * kk);
};



vector<int> Grid::ListNearbyHalos(float *X, float R)
{
	vector<int> haloIndex;

	float xMax[3], xMin[3];
	int nHalos, thisNode, indHalo, *ixMax, *ixMin;

	for (int iX = 0; iX < 3; iX++)
	{
		xMax[iX] = (X[iX] + R);
		xMin[iX] = (X[iX] - R);
		//cout << iX << "] max= " << xMax[iX] << ", min= " << xMin[iX] << endl;
	}

	ixMax = GridCoord(xMax);	
	ixMin = GridCoord(xMin);	

	for (int iX = ixMin[0]; iX <= ixMax[0]; iX++)
		for (int iY = ixMin[1]; iY <= ixMax[1]; iY++)
			for (int iZ = ixMin[2]; iZ <= ixMax[2]; iZ++)
			{
				thisNode = Index(iX, iY, iZ);		
				nHalos = haloOnGridNode[thisNode].size();
			
				for (int iH = 0; iH < nHalos; iH++)
				{
					indHalo = haloOnGridNode[thisNode][iH];

					/* Beware: indHalo can be positive (halos on the local grid) or negative (halos on the buffer) */
					haloIndex.push_back(indHalo);
				
					//if (indHalo < 0) 
					//	cout << iH << ", IndHalo: " << indHalo << endl;
				}
			}

	return haloIndex;
};



int * Grid::GridCoord(float *X)
{
	int * iX; iX = new int[3]; 
	
	for (int i = 0; i<3; i++)
		iX[i] = floor(X[i] / cellSize);

	return iX;
};



void Grid::FindNearbyNodes(int index, int nCells)
{
	int thisIndex, thisTask;
	int *iX, *jX; jX = new int[3]; 

	iX = Index2Grid(index);

	/* Periodic boundary conditions are taken care of by the Index() function, no need to implement them here */
	for (int i = -nCells; i < nCells+1; i++)
	{
		jX[0] = iX[0] + i;

		for (int j = -nCells; j < nCells+1; j++)
		{
			jX[1] = (iX[1] + j);

			for (int k = -nCells; k < nCells+1; k++)
			{ 
				jX[2] = (iX[2] + k);
				thisIndex = Index(jX[0], jX[1], jX[2]);

				/* There might be nodes shared among several tasks, so we have to loop here */
				for (int l = 0; l < globalTaskOnGridNode[thisIndex].size(); l++)
				{
						thisTask = globalTaskOnGridNode[thisIndex][l] -1;

					// Only add the task to the communication buffer if the node is not already there
					if (thisTask != locTask && thisTask > -1)
							buffNodes[thisTask].push_back(thisIndex);	
				}
			}
		}
	}
};


void Grid::SortLocNodes()
{
	sort(locNodes.begin(), locNodes.end());
	locNodes.erase(unique(locNodes.begin(), locNodes.end()), locNodes.end());
};


void Grid::AssignToGrid(float *X, int index)
{
	int *iX, thisNode = 0;
	iX = new int[3];	

	iX = GridCoord(X);
	thisNode = Index(iX[0], iX[1], iX[2]);
	haloOnGridNode[thisNode].push_back(index);	// index is positive for locHalos and negative for locBuffHalos
	taskOnGridNode[thisNode] = locTask + 1;		// Add one to distinguish from empty node (locTask = 0 has to be one, so that taskOnGridNode[] = 0 means no task is readin that)

	locNodes.push_back(thisNode);
	free(iX);
};


/* 
 * This function identifies the nodes of the buffer region for each task.
 * For each halo it determines the nodes allocated to other tasks and stores them into a list.
 * It loops on all the nodes already allocated and identifies all the neighbouring nodes within 
 * a radius of maxBufferThick size.
 * useNodes is a collection of nodes at a different snapshot. The grid might be allocated differently 
 * among the tasks at each snapshot so we compare e.g. the local nodes at 0 with those at 1 to find the buffer
 */
void Grid::FindBufferNodes(vector<int> useNodes)
{
	int nCells = 1;		// TODO find a way to select the number of buffer cells!!!!!	FIXME TODO 
				// Edit: test with 1,2 or 3 seem to give the same results...

	if (buffNodes.size() == 0)
		buffNodes.resize(totTask);

	// Do a loop on all the nodes contained in this task to find out which nodes need to be communicated
	for (int i = 0; i < useNodes.size(); i++)
		FindNearbyNodes(useNodes[i], nCells);		

	/*	Clean Buffer Nodes	*/
	for (int i = 0; i < totTask; i++)
	{
		sort(buffNodes[i].begin(), buffNodes[i].end());
		buffNodes[i].erase(unique(buffNodes[i].begin(), buffNodes[i].end()), buffNodes[i].end());
	}
};



void Grid::Info()
{

	cout << "Task=" << locTask << " has " << locNodes.size() << " loc nodes, out of " << globalTaskOnGridNode.size() << endl;	
//	int *ixMax, *ixMin;
//	ixMax = new int[3];	ixMin = new int[3];	
/*
	for (int i = ixMin[0]; i < ixMax[0]; i++)
		for (int j = ixMin[1]; j < ixMax[1]; j++)
			for (int k = ixMin[2]; k < ixMax[2]; k++)
				cout << Index(i, j, k) << " " << taskOnGridNode[Index(i, j, k)] << endl; 
*/
};


