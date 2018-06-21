#include <vector>
#include <iostream>

#include "general.h"
#include "Grid.h"

using namespace std;


Grid::Grid(int n, float size)
{
	N = n; boxSize = size;
	cellSize = (boxSize/N);

	if (locTask == 0)
		cout << "Initialized grid with N=" << N << " nodes in a box of " << boxSize << " Mpc/h" << endl; 
	
	/* This NxNxN matrix tracks the */
	//taskPerGridNode = 

};


Grid::~Grid()
{
};


