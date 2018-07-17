#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

#include "global_vars.h"
#include "utils.h"


/* 
 * 	General functions
 */
void InitLocVariables(void)
{
	// FIXME
        locVmax = 0.0; 
	totVmax = 0.0;
        locHalosSize[0] = 0; locHalosSize[1] = 0;
	nLocParts[0] = 0; nLocParts[1] = 0;

	locParts.resize(2);
	locHalos.resize(2);

	GlobalGrid[0].Init(nGrid, boxSize);
	GlobalGrid[1].Init(nGrid, boxSize);
	BufferGrid.Init(nGrid, boxSize);

	sizeHalo = sizeof(Halo);
	sizePart = sizeof(unsigned long long int);
};


unsigned int NumLines(const char * fileName)
{
	unsigned int nLines = 0;
	const char * charHead = "#";

	string lineIn;

	ifstream fileIn(fileName);

	if (fileIn)
	{
		while (getline(fileIn, lineIn))
		{
			const char * lineTmp = lineIn.c_str();

			if (lineTmp[0] != charHead[0])
				nLines++;
		}

		return nLines;
	}
	else
	{ 
		cout << "File " <<  fileName << " not found on task: " << locTask << endl;
		return -1;
	}
};


void CleanMemory(int iCat)
{

	if (locTask ==0)
		cout << "Cleaning memory for catalog " << iCat << endl;		

	locHalos[iCat].clear();
	locHalos[iCat].shrink_to_fit();

		for (int iH = 0; iH < nLocHalos[iCat]; iH++)
		{
			for (int iT = 0; iT < 6; iT++)
			{
				locParts[iCat][iH][iT].clear();
				locParts[iCat][iH][iT].shrink_to_fit();
			}
				
			locParts[iCat][iH].clear();
			locParts[iCat][iH].shrink_to_fit();
		}

		locParts[iCat].clear();
		locParts[iCat].shrink_to_fit();
		
		GlobalGrid[iCat].Clean();
};


void ShiftHalosPartsGrids()
{
	CleanMemory(0);

	if (locTask == 0)
		cout << "Shifting halos, particles and grid from 1 to 0..." << endl;

	nLocHalos[0] = nLocHalos[1];
	nLocParts[0] = nLocParts[1];

	locHalos[0] = locHalos[1];
	locParts[0].resize(nLocHalos[0]);

	for (int iH = 0; iH < nLocHalos[0]; iH++)
	{
		locParts[0][iH].resize(nPTypes);

		for (int iT = 0; iT < 6; iT++)
			locParts[0][iH][iT] = locParts[1][iH][iT];
	}
	
	GlobalGrid[0].globalTaskOnGridNode = GlobalGrid[1].globalTaskOnGridNode;
	GlobalGrid[0].taskOnGridNode = GlobalGrid[1].taskOnGridNode;
	GlobalGrid[0].locNodes = GlobalGrid[1].locNodes;
	GlobalGrid[0].buffNodes = GlobalGrid[1].buffNodes;
	GlobalGrid[0].haloOnGridNode = GlobalGrid[1].haloOnGridNode;
	GlobalGrid[0].buffOnGridNode = GlobalGrid[1].buffOnGridNode;
};


float VectorModule(float *V)
{
	return sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
};


		
// Given two halos, decide whether to compare their particle content or not
bool CompareHalos(int iHalo, int jHalo, int iOne, int iTwo)
{
	float rMax = 0.0;

	// do some check - if jHalo > nLocHalos ---> go look into the buffer halos FIXME

	rMax = locHalos[iOne][iHalo].rVir + locHalos[iTwo][jHalo].rVir; rMax *= dMaxFactor;

	// Only check for pairwise distance
	if (locHalos[iOne][iHalo].Distance(locHalos[iTwo][jHalo].X) < rMax)
		return true;
	else 
		return false;
	
};


void FindProgenitors(int iOne, int iTwo)
{
	int nStepsCounter = floor(nLocHalos[iUseCat] / 50.);
	vector<int> nCommon;

		for (int i = 0; i < nLocHalos[iOne]; i++)
		{
			if (i == nStepsCounter * floor(i / nStepsCounter) && locTask == 0)
					cout << "." << flush; 

			// FIXME instead of looping on all haloes ONLY select those within neighbouring nodes
			for (int j = 0; j < nLocHalos[iTwo]; j++)
				if (CompareHalos(i, j, iOne, iTwo))
					nCommon = CommonParticles(locParts[iOne][i], locParts[iTwo][j]);
		}
};


vector<int> CommonParticles(vector<vector<unsigned long long int>> partsHaloOne, 
	vector<vector<unsigned long long int>> partsHaloTwo)
{
	vector<int> nCommon; 
	vector<unsigned long long int>::iterator iter;
	vector<unsigned long long int> thisCommon;

	nCommon.resize(nPTypes);	

	for (int iT = 0; iT < nPTypes; iT++)
	{
		int thisSize = partsHaloOne[iT].size();
		
		if (thisSize > 0)
		{
			// This is the maximum possible number of common particles
			thisCommon.resize(thisSize);

			//cout << "Task=" << locTask << " type=" << iT << " nPart=" << partsHaloOne[iT].size() << endl;
			//cout << "Task=" << locTask << " type=" << iT << " nPart=" << partsHaloTwo[iT].size() << endl;
	
			// Find out how many particles are shared among the two arrays
			iter = set_intersection(partsHaloOne[iT].begin(), partsHaloOne[iT].end(), 
				partsHaloTwo[iT].begin(), partsHaloTwo[iT].end(), thisCommon.begin());	

			// Resize the array and free some memory
			thisCommon.resize(iter - thisCommon.begin());
			//thisCommon.shrink_to_fit();		

			// Now compute how many particles in common are there
			nCommon[iT] = thisCommon.size();
		
			//cout << "Task=" << locTask << " type=" << iT << " nCommon=" << nCommon[iT] << endl;
	
			// Clear the vector and free all the allocated memory
			thisCommon.clear();
	 		thisCommon.shrink_to_fit();
		}
	}
	
	return nCommon;

};

