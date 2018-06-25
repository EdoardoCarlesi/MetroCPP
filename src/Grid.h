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

	void RecvBufferNodes(void); 	// This function finds the nodes (on the receiving end) that need to be communicated
	void FindPatchOnTask(void);	// This uses the locXax, locXmin to determine the nodes belonging to each task	
	void AssignToGrid(float *, int);

	int Index(int, int, int);	// Given i, j, k determine their position in the array
	int * GridCoord(float *);	// Given x, y, z determine their position in the array in grid coordinates

	// Some tasks may share parts of the same node, so take care of them separately in a vector of vectors
	vector<vector<int>> globalTaskOnGridNode;	// This variable tracks the task number that is storing halos on a given grid node
	vector<int> taskOnGridNode;	// This variable tracks the task number that is storing halos on a given grid node
	vector<vector<int>> haloOnGridNode;	// For each grid node, store the (local) list of halos (indexes of locHalos)
	vector<vector<int>> buffOnGridNode;	// For each grid node, store the list of halos in the buffer (indexes of buffHalos)

};


#endif
