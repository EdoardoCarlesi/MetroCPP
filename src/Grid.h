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


#ifndef GRID_H
#define GRID_H

#include <math.h>
#include <vector>

using namespace std;


class Grid
{

public:
	Grid();
	~Grid();

	int N;
	int nNodes;

	float boxSize;	
	float cellSize;
	
	void Info(void);
	void Init(int, float);
	void Clean(void);
	
	// Sort and clean the nodes number assigned locally to each task
	void SortLocNodes(void);

	// This function determines the nodes placed on different tasks required for the buffer region
	void FindBufferNodes(vector<int>);	

	// This assigns a coordinate to the grid, and stores the id associated to the point in the node
	void AssignToGrid(float *, int);	

	// Given a node index find all the neighbouring nodes
	void FindNearbyNodes(int, int); 	

	int Index(int, int, int);	// Given i, j, k determine their position in the array
	int * Index2Grid(int);		// Given and index, i, j, k it position in grid coordinates
	int * GridCoord(float *);	// Given x, y, z determine their position in the array in grid coordinates

	/* This functions returns a list of the haloes contained within a given volume around a point inside the box */
	vector<int> ListNearbyHalos(float *, float);

	// Some tasks may share parts of the same node, so take care of them separately in a vector of vectors
	// This variable tracks the task number that is storing halos on a given grid node
	vector<vector<int>> globalTaskOnGridNode;	

	// This variable tracks the task number that is storing halos on a given grid node
	vector<int> taskOnGridNode;	

	// This variable keeps track of the grid nodes stored on the local task - it only stores the index
	vector<int> locNodes;

	// Contains a list of all tasks, for each task saves a vector containing all the indexes	
	vector<vector<int>> buffNodes;	

	// For each grid node, store the (local) list of halos (indexes of locHalos)
	vector<vector<int>> haloOnGridNode;	

	// For each grid node, store the list of halos in the buffer (indexes of buffHalos)
	vector<vector<int>> buffOnGridNode;	
};

#endif
