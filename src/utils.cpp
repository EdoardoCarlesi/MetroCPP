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
 * utils.cpp
 * This file contains simple, general pourpose functions, including most of those which are
 * responsible for cleaning up and allocating the memory on each task.
 */

#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cctype>
#include <map>

#include "global_vars.h"
#include "utils.h"

using namespace std;


/* 
 * 	General functions
 */
void InitLocVariables(void)
{
        locVmax = 0.0; 
	totVmax = 0.0;
        locHalosSize[0] = 0; locHalosSize[1] = 0;
	nLocParts[0] = 0; nLocParts[1] = 0;

	nPTypes = NPTYPES;

	locParts.resize(2);
	locHalos.resize(2);
	locMTrees.resize(2);

	locMapParts.resize(2);

#ifndef ZOOM
	GlobalGrid[0].Init(nGrid, boxSize);
	GlobalGrid[1].Init(nGrid, boxSize);
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

	if (locHalos[iCat].size() > 0)
	{
		locHalos[iCat].clear();
		locHalos[iCat].shrink_to_fit();
	}

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

	nextMapTrees.clear();
	thisMapTrees.clear();
	locMapParts[iCat].clear();

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

	//cout << locTask << " Orphan: " << locOrphHalos.size() << endl;

#ifndef ZOOM
	/* Keep track of the orphan halos at the next step */
	for (int iO = 0; iO < locOrphHalos.size(); iO++)
	{
		//if (!locOrphHalos[iO].isToken || locOrphHalos[iO].nOrphanSteps == 0)	// FIXME nOrphanSteps == 0 first, why?
		//	cout << "Bad assignment of token status to halo " << iO << " on task " << locTask << endl;

		if (!locOrphHalos[iO].isToken) 
			locOrphHalos[iO].isToken = true;

		if (locOrphHalos[iO].nOrphanSteps == 0)	
			locOrphHalos[iO].isToken += 1;

		locHalos[0].push_back(locOrphHalos[iO]);
	}
#endif

	if (runMode == 0 || runMode == 2)
	{ 
		locParts[0].resize(nLocHalos[0]);
		nLocParts[0] = nLocParts[1];

		for (int iH = 0; iH < nLocHalos[0]; iH++)
		{
			locParts[0][iH].resize(nPTypes);

			for (int iT = 0; iT < nPTypes; iT++)
			{
				locParts[0][iH][iT].swap(locParts[1][iH][iT]);

				for (auto const partID : locParts[0][iH][iT])
				{
					Particle thisParticle;
 	                	        thisParticle.haloID = locHalos[0][iH].ID;
                                	thisParticle.type   = iT; 
                                	locMapParts[0][partID].push_back(thisParticle);
				}

			}
		}

#ifdef ZOOM
	}	
#else
		for (int iO = 0; iO < locOrphHalos.size(); iO++)
		{
			locParts[0].push_back(locOrphParts[iO]);
			locParts[0][nLocHalos[0]+iO].resize(nPTypes);

			for (int iT = 0; iT < nPTypes; iT++)
			{
				locParts[0][nLocHalos[0]+iO][iT].swap(locOrphParts[iO][iT]);
				for (auto const & partID : locParts[0][nLocHalos[0]+iO][iT])
				{
					Particle thisParticle;
 	                	        thisParticle.haloID = locOrphHalos[iO].ID;

                                	thisParticle.type   = iT;
                                	locMapParts[0][partID].push_back(thisParticle);
				}
			}
		}	
	}	// runMode 0 or 2
		
	/* Now free and reset the orphan halo trackers */
	nLocHalos[0] += locOrphHalos.size();
	locOrphHalos.clear();
	locOrphHalos.shrink_to_fit();

	/* Re assign the halos to the locNodes on GlobalGrid[0] 
	 * DO NOT COPY IT FROM GlobalGrid[1] - this contains also the buffer nodes in locNodes, it is difficult 
	 * to disentangle, and will grow the buffer exponentially at each step */
	
	GlobalGrid[0].Init(nGrid, boxSize);

	for (int iH = 0; iH < nLocHalos[0]; iH++)
		GlobalGrid[0].AssignToGrid(locHalos[0][iH].X, iH);

	/* Now clean the halo & particle buffers */
	locBuffHalos.clear();
	locBuffHalos.shrink_to_fit();
	
	if (runMode == 0 || runMode == 2)
	{

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
	}
	
#endif	// ZOOM mode
	
	if (locTask == 0)
		cout << "Grid, halo and particle data has been cleaned and copied 1 ---> 0." << endl;

	CleanMemory(1);

#ifndef ZOOM
	/* Reallocate a grid for the next loop */
	GlobalGrid[1].Init(nGrid, boxSize);
#endif
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
	vector<float> newVec;
	map<float, int> floatPos;
	idx.resize(nVec);

	newVec = vec;

	for (int iV = 0; iV < nVec; iV++)
		floatPos.insert(make_pair(vec[iV], iV));

	sort(newVec.begin(), newVec.end());

	float thisG = vec[0];

	for (int iV = 0; iV < nVec; iV++)
		idx[iV] = floatPos[newVec[iV]];

	/*
	
	for (int iV = 0; iV < nVec; iV++)
		idx[iV] = iV;

	sort(idx.begin(), idx.end(), [&vec](int i1, int i2) 
		{return vec[i1] < vec[i2];});

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

