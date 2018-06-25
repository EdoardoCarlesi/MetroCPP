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
	float boxSize;	
	float cellSize;
	
	void Info(void);
	void Init(int, float);

	void RecvBufferNodes(void); 	// This function finds the nodes (on the receiving end) that need to be communicated
	void FindPatchOnTask(void);	// This uses the locXax, locXmin to determine the nodes belonging to each task	

	int Index(int, int, int);	// Given i, j, k determine their position in the array
	int * GridCoord(float *);	// Given x, y, z determine their position in the array in grid coordinates

	vector<int> taskOnGridNode;

};


#endif
