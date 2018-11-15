#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cctype>

#include "global_vars.h"
#include "utils.h"


/* 
 * 	General functions
 */
void InitLocVariables(void)
{
        locVmax = 0.0; 
	totVmax = 0.0;
        locHalosSize[0] = 0; locHalosSize[1] = 0;
	nLocParts[0] = 0; nLocParts[1] = 0;

	minPartCmp = 10;

	nPTypes = NPTYPES;

	locParts.resize(2);
	locHalos.resize(2);
	locMTrees.resize(2);

#ifndef ZOOM
	GlobalGrid[0].Init(nGrid, boxSize);
	GlobalGrid[1].Init(nGrid, boxSize);
	BufferGrid.Init(nGrid, boxSize);
#endif

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
		exit(0);
	}
};



void CleanMemory(int iCat)
{

	if (locTask ==0)
		cout << "Cleaning memory for catalog " << iCat << endl;	

	locHalos[iCat].clear();
	locHalos[iCat].shrink_to_fit();
	nLocHalos[iCat] = 0;

	/* Clean the particles if not running in post processing mode only */
	if (runMode == 0 || runMode == 2)
	{
		for (int iH = 0; iH < nLocHalos[iCat]; iH++)
		{
			for (int iT = 0; iT < nPTypes; iT++)
			{
				locParts[iCat][iH][iT].clear();
				locParts[iCat][iH][iT].shrink_to_fit();
			}
				
			locParts[iCat][iH].clear();
			locParts[iCat][iH].shrink_to_fit();
		}

		locParts[iCat].clear();
		locParts[iCat].shrink_to_fit();
	}

#ifndef ZOOM		
	GlobalGrid[iCat].Clean();
#endif
};



vector<string> SplitString (string strIn, string delim)
{
	string whiteSpace, tabSpace, endSpace; 
	vector<string> results, cleanResults;
	int cutAt, iChar = 0;

	whiteSpace = " "; tabSpace = '\t'; endSpace = '\n'; 

	while((cutAt = strIn.find_first_of(delim)) != strIn.npos )
	{
		if(cutAt > 0)
        		results.push_back(strIn.substr(0, cutAt));
	
    		strIn = strIn.substr(cutAt+1);
    	}

	if(strIn.length() > 0)
		results.push_back(strIn);

	for (int iR = 0; iR < results.size(); iR++)
	{
		string newStr = results[iR];
		newStr.erase(remove(newStr.begin(), newStr.end(), ' '), newStr.end());
		newStr.erase(remove(newStr.begin(), newStr.end(), '\t'), newStr.end());
		newStr.erase(remove(newStr.begin(), newStr.end(), " "), newStr.end());
		cleanResults.push_back(newStr);
	}

    return cleanResults;
}



void ShiftHalosPartsGrids()
{
	CleanMemory(0);

	if (locTask == 0)
		cout << "Shifting halos, particles and grid from 1 to 0..." << endl;

	nLocHalos[0] = nLocHalos[1];
	locHalos[0] = locHalos[1];

	if (runMode == 0 || runMode == 2)
	{ 
		locParts[0].resize(nLocHalos[0]);
		nLocParts[0] = nLocParts[1];
		//cout << nLocHalos[1] << ", " << locParts[0].size() << ", " << locParts[1].size() << endl;

		for (int iH = 0; iH < nLocHalos[0]; iH++)
		{
			locParts[0][iH].resize(nPTypes);

			for (int iT = 0; iT < nPTypes; iT++)
				locParts[0][iH][iT].swap(locParts[1][iH][iT]);
		}
	}
	
#ifndef ZOOM
	GlobalGrid[0].taskOnGridNode.swap(GlobalGrid[1].taskOnGridNode);
	GlobalGrid[0].locNodes.swap(GlobalGrid[1].locNodes);

	int sizeTaskOnGridNode = GlobalGrid[1].globalTaskOnGridNode.size();
	GlobalGrid[0].globalTaskOnGridNode.resize(sizeTaskOnGridNode);
	for (int iN = 0; iN < sizeTaskOnGridNode; iN++)
		GlobalGrid[0].globalTaskOnGridNode[iN].swap(GlobalGrid[1].globalTaskOnGridNode[iN]);

	int sizeBuffNodes = GlobalGrid[1].buffNodes.size();
	GlobalGrid[0].buffNodes.resize(sizeBuffNodes);
	for (int iN = 0; iN < sizeBuffNodes; iN++)
		GlobalGrid[0].buffNodes[iN].swap(GlobalGrid[1].buffNodes[iN]);

	int sizeHaloOnGridNode = GlobalGrid[1].haloOnGridNode.size();
	GlobalGrid[0].haloOnGridNode.resize(sizeHaloOnGridNode);
	for (int iN = 0; iN < sizeHaloOnGridNode; iN++)
		GlobalGrid[0].haloOnGridNode[iN].swap(GlobalGrid[1].haloOnGridNode[iN]);

	int sizeBuffOnGridNode = GlobalGrid[1].buffOnGridNode.size();
	GlobalGrid[0].buffOnGridNode.resize(sizeBuffOnGridNode);
	for (int iN = 0; iN < sizeBuffOnGridNode; iN++)
		GlobalGrid[0].buffOnGridNode[iN].swap(GlobalGrid[1].buffOnGridNode[iN]);
#endif
	
	/* Now clean the halo & particle buffers */
	locBuffHalos.clear();
	locBuffHalos.shrink_to_fit();
	
	for (int iP = 0; iP < locBuffParts.size(); iP++)
	{
		for (int iT = 0; iT < nPTypes; iT++)
		{
			locBuffParts[iP][iT].clear();
			locBuffParts[iP][iT].shrink_to_fit();
		}

		locBuffParts[iP].clear();
		locBuffParts[iP].shrink_to_fit();
	}
	
	locBuffParts.clear();
	locBuffParts.shrink_to_fit();

	/* Data from sim 1 has been copied into 0, so free 1 */
	CleanMemory(1);
	
	if (locTask == 0)
		cout << "Grid, halo and particle data has been copied and cleaned." << endl;

	/* Reallocate a grid for the next loop */
	GlobalGrid[1].Init(nGrid, boxSize);
};



float VectorModule(float *V)
{
	return sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
};



float *UnitVector(float *V)
{
	float *nV; nV = (float *) calloc(3, sizeof(float));
	float normV; normV = VectorModule(V);

	for (int iN = 0; iN < 3; iN++)
	{
		nV[iN] = V[iN] / normV;
	}

	return nV;
};



//int SortIndexes(vector<float> vec) {
vector<int> SortIndexes(vector<float> vec) {
	int nVec = vec.size();

	vector<int> idx;
	idx.resize(nVec);
	
	for (int iV = 0; iV < nVec; iV++)
		idx[iV] = iV;

	sort(idx.begin(), idx.end(), [&vec](int i1, int i2) 
		{return vec[i1] > vec[i2];});

	/*
	//int idx0 = 0; 
	//float maxV = 100000.0;
	//float maxV = 0.0;
	//for (int iD = 0; iD < nVec; iD++)
	//	cout << "* " << iD << " vec=" << vec[iD] << " idx=" << idx[iD] << endl;
	//sort(idx.begin(), idx.end(), [&vec](int i1, int i2)
	//	{return vec[i1] < vec[i2];});

	for (int iD = 0; iD < nVec; iD++)
	{
		if (vec[iD] > maxV)
		{
			maxV = vec[iD];
			idx0 = iD;
		}
	}
	//return idx0;
	*/

	return idx;
};

