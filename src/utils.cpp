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
	minPartHalo = 30;

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
		return -1;
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

		for (int iH = 0; iH < nLocHalos[0]; iH++)
		{
			locParts[0][iH].resize(nPTypes);

			for (int iT = 0; iT < nPTypes; iT++)
				locParts[0][iH][iT] = locParts[1][iH][iT];
		}
	}	

#ifndef ZOOM
	GlobalGrid[0].globalTaskOnGridNode = GlobalGrid[1].globalTaskOnGridNode;
	GlobalGrid[0].taskOnGridNode = GlobalGrid[1].taskOnGridNode;
	GlobalGrid[0].locNodes = GlobalGrid[1].locNodes;
	GlobalGrid[0].buffNodes = GlobalGrid[1].buffNodes;
	GlobalGrid[0].haloOnGridNode = GlobalGrid[1].haloOnGridNode;
	GlobalGrid[0].buffOnGridNode = GlobalGrid[1].buffOnGridNode;
#endif
	
	/* Data from sim 1 has been copied into 0, so free 1 */
	CleanMemory(1);

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

