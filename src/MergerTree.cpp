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
 * MergerTree.cpp:
 * This file holds the core routines and defines the classes that are required to consistently compute 
 * the merger trees.
 */

#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <map>

#include "MergerTree.h"
#include "Halo.h"
#include "utils.h"
#include "global_vars.h"

using namespace std;



HaloTree::HaloTree()
{
};



HaloTree::~HaloTree()
{
	Clean();
};



void HaloTree::Clean()
{
	for (int iM = 0; iM < progHalo.size(); iM++)
			progHalo[iM].clear();
	
	mainHalo.clear();
	progHalo.clear();
};



void MergerTree::Info()
{
	cout << "Task=" << locTask << " " << mainHalo.ID << " " << idProgenitor.size() << endl;

	if (isOrphan)
		for (int iP = 0; iP < idProgenitor.size(); iP++)
			cout << "Ind=" << indexProgenitor[iP] << " ID= " << idProgenitor[iP] << endl;
	else
		for (int iP = 0; iP < idProgenitor.size(); iP++)
			cout << "Ind=" << indexProgenitor[iP] << " ID= " << idProgenitor[iP] << " NP=" << nCommon[1][iP] << endl;

};



void MergerTree::Clean()
{
	for (int iC = 0; iC < nCommon.size(); iC++)
		nCommon[iC].clear();

	nCommon.clear();

	idProgenitor.clear();
	indexProgenitor.clear();
	progHalos.clear();
};



MergerTree::MergerTree()
{
};

MergerTree::~MergerTree()
{
};


void MergerTree::AssignMap()
{
	int nProgs = 0, nCommTot = 0, iP = 0;
	map<uint64_t, vector<int>>::iterator thisMap;
	
	nCommon.resize(nPTypes);
	
	/* Each merger tree stores the halo ids in a map that stores the number of particles shared with each progenitor */
	for (thisMap = indexCommon.begin(); thisMap != indexCommon.end(); thisMap++)
	{
		uint64_t thisHaloID = thisMap->first;
		vector<int> thisCommon = thisMap->second;

		nCommTot = 0;

		for (int iT = 0; iT < nPTypes; iT++)
			nCommTot += thisCommon[iT];

		if (nCommTot > minPartCmp)	
		{
			for (int iT = 0; iT < nPTypes; iT++)
				nCommon[iT].push_back(thisCommon[iT]);
			
			idProgenitor.push_back(thisHaloID);
			nProgs++;
		}

		iP++;
	}

	if (nProgs == 0)
	{
		isOrphan = true;

	} else {
		/* We allocate these vectors here, the indexes and the halos will be copied inside later on */
		isOrphan = false;
		indexProgenitor.resize(nProgs);
		progHalos.resize(nProgs);
	}
};


/* This sorting algorithm might be very inefficient but it's straightforward to implement, 
 * plus we will rarely deal with halos with more than 10^3 progenitors */
void MergerTree::SortByMerit()
{
	vector<uint64_t> tmpIdx;
	vector<vector<int>> tmpNCommon;
	vector<float> allMerit;
	vector<int> idx, tmpIndex;
	vector<Halo> tmpProgHalos;
	float merit = 0.0;

	for (int iM = 0; iM < idProgenitor.size(); iM++)
	{
		int nComm = 0;
		double ratioM = 0;

		for (int iC = 0; iC < nPTypes; iC++)
			nComm += nCommon[iC][iM];
		
		int nPH0 = 0, nPH1 = 0;
 
		for (int iP = 0; iP < nPTypes; iP++)
		{	
			nPH0 += mainHalo.nPart[iP];
			nPH1 += progHalos[iM].nPart[iP];
		}

		//ratioM = (float) mainHalo.nPart[1] / (float) progHalos[iM].nPart[1];
		ratioM = (float) nPH0 / (float) nPH1;
	
		if (ratioM < 1.0) ratioM = 1.0 / ratioM;

		merit = nComm  / (ratioM*1.0001 - 1.0);
		merit *= (1.0 + 0.00001 * iM);	// We change the merit slightly, 
						// to avoid confusion when two halos have the same number of particles 
						// and the same number of particles shared with the host halo

		allMerit.push_back(merit);
	}
	
	idx = SortIndexes(allMerit);
	tmpProgHalos.resize(idx.size());
	tmpIndex.resize(idx.size());
	tmpIdx.resize(idx.size());
	tmpNCommon.resize(nPTypes);

	for (int iT = 0; iT < nPTypes; iT++)
		tmpNCommon[iT].resize(idx.size());

	for (int iM = 0; iM < idProgenitor.size(); iM++)
	{
		/* Inverse sort - from largest to smallest */
		int oldIdx, newIdx;
		oldIdx = idProgenitor.size() - iM - 1;
		newIdx = idx[iM];

		for (int iT = 0; iT < nPTypes; iT++)
			tmpNCommon[iT][oldIdx] = nCommon[iT][newIdx];

		tmpIdx[oldIdx] = idProgenitor[newIdx];
		tmpIndex[oldIdx] = indexProgenitor[newIdx];
		tmpProgHalos[oldIdx] = progHalos[newIdx];
	}	

	for (int iM = 0; iM < idProgenitor.size(); iM++)
	{
		idProgenitor[iM] = tmpIdx[iM];
		indexProgenitor[iM] = tmpIndex[iM];
		progHalos[iM] = tmpProgHalos[iM];

		for (int iT = 0; iT < nPTypes; iT++)
			nCommon[iT][iM] = tmpNCommon[iT][iM];
	}

	tmpIdx.clear();
	tmpIdx.shrink_to_fit();
	tmpIndex.clear();
	tmpIndex.shrink_to_fit();
	tmpProgHalos.clear();
	tmpProgHalos.shrink_to_fit();
};

	       /****************************************************************************
		* These are general functions which do not belong to the class Merger Tree *
		****************************************************************************/

#ifndef ZOOM		

/* This is a very fast way of comparing particles content of halos across snapshots, that relies on maps 
 * and scales linearly with the number of particles on each task */

void FindProgenitors(int iOne, int iTwo)
{
	int nLoopHalos[2], iOldOrphans = 0, iFixOrphans = 0, nLocOrphans = 0, nLocTrees = 0; 

	/* Loop also on the buffer halos, in the backward loop only! */
	if (iOne == 1)
	{
		nLoopHalos[iOne] = nLocHalos[iOne] + locBuffHalos.size();
		nLoopHalos[iTwo] = nLocHalos[iTwo]; 
	} else {
		nLoopHalos[iOne] = nLocHalos[iOne];
		nLoopHalos[iTwo] = nLocHalos[iTwo] + locBuffHalos.size();
	}

	locMTrees[iOne].clear();
	locMTrees[iOne].shrink_to_fit();
	locMTrees[iOne].resize(nLoopHalos[iOne]);

	Halo thisHalo;

#ifdef VERBOSE
	if (locTask == 0)
	{
		cout << iOne << ", Loop, " << nLocHalos[iOne] << ",  iTwo " << iTwo << " " << nLocHalos[iTwo] << endl;
		cout << iOne << ", Loc , " << nLoopHalos[iOne] << ",  iTwo " << iTwo << " " << nLoopHalos[iTwo] << endl;
	}
#endif

	/* Reset the tree maps for the inverse comparison */
	if (iOne == 1)
	{
		thisMapTrees.clear();
		nextMapTrees.clear();
	}

	/* Here we initialize (again, for safety and because it's cheap) the maps of halo indexes within the locHalos vectors
	 * and their IDs for faster identification in the loop on the particles */
	for (int iL = 0; iL < nLoopHalos[iOne]; iL++)
	{
		int iH = 0;

		if (iL < nLocHalos[iOne])
		{	
			iH = iL;
			thisHalo = locHalos[iOne][iH];
		} else {
			iH = nLocHalos[iOne]-iL-1;
			thisHalo = locBuffHalos[-iH-1];
		}

		/* Here we link every halo id to its position in the locHalo vector */
		thisMapTrees[thisHalo.ID] = iL;

		locMTrees[iOne][iL].mainHalo = thisHalo; 
		locMTrees[iOne][iL].nCommon.resize(nPTypes);
	}

	/* Map the iTwo halo IDs to their indexes */
	for (int iL = 0; iL < nLoopHalos[iTwo]; iL++)
	{
		int iH = 0;

		if (iL < nLocHalos[iTwo])
		{	
			iH = iL;
			thisHalo = locHalos[iTwo][iH];
		} else {
			iH = nLocHalos[iTwo]-iL-1;
			thisHalo = locBuffHalos[-iH-1];
		}
	
		/* This map connects halo IDs & their indexes in the SECOND locHalo structure */
		nextMapTrees[thisHalo.ID] = iL;
	}

	/* Here we loop on all the particles, each particle keeps track of the Halos it belongs to. 
	 * We match particle IDs in iOne with particle IDs in iTwo, and count the total number of 
	 * particles shared by their two host halos. */
	for (auto const& thisMap : locMapParts[iOne]) 
	{
		uint64_t thisID = thisMap.first;		// This particle ID
		vector<Particle> thisParticle = thisMap.second;	 	// How many halos (halo IDs) share this particle

		/* Loop on the halos on iOne to which this particle belongs */
		for (int iH = 0; iH < thisParticle.size(); iH++)
		{
			vector<Particle> nextParticle = locMapParts[iTwo][thisID];
			int thisTreeIndex = thisMapTrees[thisParticle[iH].haloID];
	
			/* This same particle on iTwo is also shared by some halos: do the match with the haloIDs on iOne. */
			for (int iN = 0; iN < nextParticle.size(); iN++)
			{
				uint64_t nextHaloID = nextParticle[iN].haloID;

				/* If this ID is not in the list of progenitor IDs, then initialize the indexCommon 
				   map and initialize the number of common particles */
				if (locMTrees[iOne][thisTreeIndex].indexCommon.find(nextHaloID) == 
					locMTrees[iOne][thisTreeIndex].indexCommon.end())
				{
					locMTrees[iOne][thisTreeIndex].indexCommon[nextHaloID].resize(nPTypes);
					locMTrees[iOne][thisTreeIndex].indexCommon[nextHaloID][nextParticle[iN].type] = 1;
				} else {	/* If the Halo ID is already in the index of the halos with common particles,
						   then add ++ to the particle type shared */ 
					locMTrees[iOne][thisTreeIndex].indexCommon[nextHaloID][nextParticle[iN].type]++;
				}
			} // Loop on the halos in the iTwo particles
		} // Loop on the halos in the iOne particles  
	} // Loop on all the iOne particles

	int iOrph = 0;

	/* Once the first loop on the particles is done, we need to fix ALL merger trees
	   moving all the data stored in the map to the "standard" index & id vectors */
	for (int iM = 0; iM < locMTrees[iOne].size(); iM++)
		locMTrees[iOne][iM].AssignMap();

	/* Now clean and reconstruct the local merger trees */
	for (int iM = 0; iM < locMTrees[iOne].size(); iM++)
	{
		uint64_t thisHaloID;
		int nProgs = locMTrees[iOne][iM].idProgenitor.size();
		int thisHaloIndex = 0;

		for (int iP = 0; iP < nProgs; iP++)
		{
			thisHaloID = locMTrees[iOne][iM].idProgenitor[iP];
			thisHaloIndex = nextMapTrees[thisHaloID];			

			/* Here we fill the indexProgenitor & progHalos */
			locMTrees[iOne][iM].indexProgenitor[iP] = thisHaloIndex;

			if (thisHaloIndex >= nLocHalos[iTwo])
				locMTrees[iOne][iM].progHalos[iP] = locBuffHalos[thisHaloIndex-nLocHalos[iTwo]];
			else
				locMTrees[iOne][iM].progHalos[iP] = locHalos[iTwo][thisHaloIndex];
		}

		locMTrees[iOne][iM].SortByMerit();

		/* Orphan halos are identified in the forward search only */
		if (iOne == 0)
		{

#ifdef NOPTYPE
			if (locMTrees[iOne][iM].idProgenitor.size() == 0 && 
				locMTrees[iOne][iM].mainHalo.nPart[0] > minPartHalo)
#else
			if (locMTrees[iOne][iM].idProgenitor.size() == 0 && 
				locMTrees[iOne][iM].mainHalo.nPart[1] > minPartHalo)
#endif
			{
				Halo thisHalo = locHalos[iOne][iM];
				thisHalo.isToken = true;
				thisHalo.nOrphanSteps++;

				if (thisHalo.nOrphanSteps > 1)
					iOldOrphans++;

				/* Update the container of local orphan halos */
				locOrphIndex.push_back(iM);
				locOrphHalos.push_back(thisHalo);

				/* Update the local mtree with a copy of itself */
				locMTrees[iOne][iM].isOrphan = true;
				locMTrees[iOne][iM].idProgenitor.push_back(locHalos[iOne][iM].ID);

				/* Update the particle content */
				locOrphParts.push_back(locParts[iOne][iM]);
				locOrphParts[nLocOrphans].resize(nPTypes);

				/* Copy particle blocks divided by particle type */
				for (int iP = 0; iP < nPTypes; iP++)
					copy(locParts[iOne][iM][iP].begin(), locParts[iOne][iM][iP].begin(), 
						back_inserter(locOrphParts[nLocOrphans][iP]));

				nLocOrphans++;
			} else { // This halo has a progenitor

				if (locHalos[iOne][iM].isToken)
					iFixOrphans++;

				locMTrees[iOne][iM].isOrphan = false;
				nLocTrees++;
			}
		} 
	}

	/* Trace the orphans in the forward loop */
	if (iOne == 0)
	{
		nLocOrphans = locOrphHalos.size();
		int nTotOrphans = 0, nTotFix = 0, nTotOld = 0, nTotTrees = 0; 

		//cout << "\nFound " << nLocOrphans << " orphan halos on task " << locTask << ", " << locOrphIndex.size() << endl;
		MPI_Reduce(&nLocTrees,   &nTotTrees, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&nLocOrphans, &nTotOrphans, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&iOldOrphans, &nTotOld, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&iFixOrphans, &nTotFix, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

		if (locTask == 0)
		{	
			cout << endl;
			cout << "  Tracking a total of " << nTotTrees << " halos on " << totTask << " tasks. "  << endl;
			cout << "  On all tasks, there are: " << endl;
			cout << "     ---> " << nTotOrphans << " total orphan halos." << endl;
			cout << "     ---> " << nTotOrphans - nTotOld << " new orphan halos." << endl;
			cout << "     ---> " << nTotOld << " orphans since more than one step." << endl;
			cout << "     ---> " << nTotFix << " orphans reconnected to their descendants.  " << endl;
			cout << "  On the master task there are " << nLocTrees << " halos as well as: " << endl;
			cout << "     ---> " << nLocOrphans << " total orphan halos." << endl;
			cout << "     ---> " << nLocOrphans - iOldOrphans << " new orphan halos." << endl;
			cout << "     ---> " << iOldOrphans << " orphans since more than one step." << endl;
			cout << "     ---> " << iFixOrphans << " orphans reconnected to their descendants.  " << endl;
		}
	}

};		/* End of the find progenitor function in full box mode */

#else		 /* The Find Progenitors function in ZOOM MODE is different than the standard one */

void FindProgenitors(int iOne, int iTwo)
{
	int iOldOrphans = 0, iFixOrphans = 0, nLocOrphans = 0, nLocTrees = 0; 

	locMTrees[iOne].clear();
	locMTrees[iOne].shrink_to_fit();
	locMTrees[iOne].resize(nLocHalos[iOne]);

	Halo thisHalo;

#ifdef VERBOSE
	if (locTask == 0)
		cout << iOne << ", Loop, " << nLocHalos[iOne] << ",  iTwo " << iTwo << " " << nLocHalos[iTwo] << endl;
#endif

	/* Reset the tree maps for the inverse comparison */
	if (iOne == 1)
	{
		thisMapTrees.clear();
		nextMapTrees.clear();
	}

	/* Here we initialize (again, for safety and because it's cheap) the maps of halo indexes within the locHalos vectors
	 * and their IDs for faster identification in the loop on the particles */
	for (int iL = 0; iL < nLocHalos[iOne]; iL++)
	{
		thisHalo = locHalos[iOne][iL];

		/* Here we link every halo id to its position in the locHalo vector */
		thisMapTrees[thisHalo.ID] = iL;

		locMTrees[iOne][iL].mainHalo = thisHalo; 
		locMTrees[iOne][iL].nCommon.resize(nPTypes);
	}

	/* Map the iTwo halo IDs to their indexes */
	for (int iL = 0; iL < nLocHalos[iTwo]; iL++)
	{
		thisHalo = locHalos[iTwo][iL];
	
		/* This map connects halo IDs & their indexes in the SECOND locHalo structure */
		nextMapTrees[thisHalo.ID] = iL;
	}

	/* Here we loop on all the particles, each particle keeps track of the Halos it belongs to. 
	 * We match particle IDs in iOne with particle IDs in iTwo, and count the total number of 
	 * particles shared by their two host halos. */
	for (auto const& thisMap : locMapParts[iOne]) 
	{
		uint64_t thisID = thisMap.first;		// This particle ID
		vector<Particle> thisParticle = thisMap.second;	 	// How many halos (halo IDs) share this particle

		/* Loop on the halos on iOne to which this particle belongs */
		for (int iH = 0; iH < thisParticle.size(); iH++)
		{
			vector<Particle> nextParticle = locMapParts[iTwo][thisID];
			int thisTreeIndex = thisMapTrees[thisParticle[iH].haloID];
	
			/* This same particle on iTwo is also shared by some halos: do the match with the haloIDs on iOne. */
			for (int iN = 0; iN < nextParticle.size(); iN++)
			{
				uint64_t nextHaloID = nextParticle[iN].haloID;

				/* If this ID is not in the list of progenitor IDs, then initialize the indexCommon 
				   map and initialize the number of common particles */
				if (locMTrees[iOne][thisTreeIndex].indexCommon.find(nextHaloID) == 
					locMTrees[iOne][thisTreeIndex].indexCommon.end())
				{
					locMTrees[iOne][thisTreeIndex].indexCommon[nextHaloID].resize(nPTypes);
					locMTrees[iOne][thisTreeIndex].indexCommon[nextHaloID][nextParticle[iN].type] = 1;
				} else {	/* If the Halo ID is already in the index of the halos with common particles,
						   then add ++ to the particle type shared */ 
					locMTrees[iOne][thisTreeIndex].indexCommon[nextHaloID][nextParticle[iN].type]++;
				}
			} // Loop on the halos in the iTwo particles
		} // Loop on the halos in the iOne particles  
	} // Loop on all the iOne particles

	int iOrph = 0;

	/* Once the first loop on the particles is done, we need to fix ALL merger trees
	   moving all the data stored in the map to the "standard" index & id vectors */
	for (int iM = 0; iM < locMTrees[iOne].size(); iM++)
		locMTrees[iOne][iM].AssignMap();

	/* Now clean and reconstruct the local merger trees */
	for (int iM = 0; iM < locMTrees[iOne].size(); iM++)
	{
		uint64_t thisHaloID;
		int nProgs = locMTrees[iOne][iM].idProgenitor.size();
		int thisHaloIndex = 0;

		for (int iP = 0; iP < nProgs; iP++)
		{
			thisHaloID = locMTrees[iOne][iM].idProgenitor[iP];
			thisHaloIndex = nextMapTrees[thisHaloID];			

			/* Here we fill the indexProgenitor & progHalos */
			locMTrees[iOne][iM].indexProgenitor[iP] = thisHaloIndex;
			locMTrees[iOne][iM].progHalos[iP] = locHalos[iTwo][thisHaloIndex];
		}

		locMTrees[iOne][iM].SortByMerit();

		/* Orphan halos are identified in the forward search only */
		if (iOne == 0)
		{
			if (locMTrees[iOne][iM].idProgenitor.size() == 0 && 
				locMTrees[iOne][iM].mainHalo.nPart[1] > minPartHalo)
			{
				Halo thisHalo = locHalos[iOne][iM];
				thisHalo.isToken = true;
				thisHalo.nOrphanSteps++;

				if (thisHalo.nOrphanSteps > 1)
					iOldOrphans++;

				/* Update the container of local orphan halos */
				locOrphIndex.push_back(iM);
				locOrphHalos.push_back(thisHalo);

				/* Update the local mtree with a copy of itself */
				locMTrees[iOne][iM].isOrphan = true;
				locMTrees[iOne][iM].idProgenitor.push_back(locHalos[iOne][iM].ID);

				/* Update the particle content */
				locOrphParts.push_back(locParts[iOne][iM]);
				locOrphParts[nLocOrphans].resize(nPTypes);

				/* Copy particle blocks divided by particle type */
				for (int iP = 0; iP < nPTypes; iP++)
					copy(locParts[iOne][iM][iP].begin(), locParts[iOne][iM][iP].begin(), 
						back_inserter(locOrphParts[nLocOrphans][iP]));

				nLocOrphans++;
			} else { // This halo has a progenitor

				if (locHalos[iOne][iM].isToken)
					iFixOrphans++;

				locMTrees[iOne][iM].isOrphan = false;
				nLocTrees++;
			}
		} 
	}

	/* Stats about orphan halos in the fwd loop */
	if (iOne == 0)
	{
		cout << endl;
		cout << "  On the master task there are " << nLocTrees << " halos as well as: " << endl;
		cout << "     ---> " << nLocOrphans << " total orphan halos." << endl;
		cout << "     ---> " << nLocOrphans - iOldOrphans << " new orphan halos." << endl;
		cout << "     ---> " << iOldOrphans << " orphans since more than one step." << endl;
		cout << "     ---> " << iFixOrphans << " orphans reconnected to their descendants.  " << endl;
	}

};
#endif		// ifndef ZOOM


void InitTrees(int nUseCat)
{
	locCleanTrees.resize(nUseCat-1);
};


/* This function compares the forward/backward connections to determine the unique descendant of each halo */
void CleanTrees(int iStep)
{
	int thisIndex = 0, halosPerTask = 0, halosRemaind = 0;
	int nErr = 0;
        halosPerTask = int(nLocHalos[0] / totTask);
        halosRemaind = nLocHalos[0] % totTask;

#ifdef ZOOM
	if (locTask < halosRemaind)
		halosPerTask += 1;
#else
	halosPerTask = nLocHalos[0]; 
#endif

	if (locTask == 0)
		cout << "Cleaning Merger Tree connections for " << halosPerTask << " halos." << endl;

#ifdef VERBOSE
	if (locTask == 0)
		cout << "nHalos " << locHalos[0].size() << ", nTrees: " << locMTrees[0].size() << endl; 
#endif

	for (int kTree = 0; kTree < halosPerTask; kTree++)
	{

#ifdef ZOOM	// In ZOOM mode each task holds all the halos, but only analyzes some of them
		int iTree = locTask + kTree * totTask;
#else
		int iTree = kTree;
#endif
		uint64_t mainID = locHalos[0][iTree].ID;
		int nProgSize = locMTrees[0][iTree].idProgenitor.size();

		MergerTree mergerTree;
		mergerTree.nCommon.resize(nPTypes);
		mergerTree.mainHalo = locHalos[0][iTree];
		mergerTree.isOrphan = locMTrees[0][iTree].isOrphan;

		if (mergerTree.isOrphan)
		{
			nProgSize = 0;
			locMTrees[0][iTree].idProgenitor[0] = mainID;
			mergerTree.idProgenitor.push_back(mainID);
			mergerTree.progHalos.push_back(locHalos[0][iTree]);
			
			for (int iT = 0; iT < nPTypes; iT++)
			{
				mergerTree.nCommon[iT].resize(1);
				mergerTree.nCommon[iT][0] = locHalos[0][iTree].nPart[iT];
			}
		}

		/* At each step we only record the connections between halos in catalog 0 and catalog 1, without attempting at a
		 * reconstruction of the full merger history. This will be done later. */
		for (int iProg = 0; iProg < nProgSize; iProg++)
		{
			Halo progHalo;	
			int jTree = locMTrees[0][iTree].indexProgenitor[iProg];
			uint64_t progID = locMTrees[0][iTree].idProgenitor[iProg];
			uint64_t descID;

#ifndef ZOOM 		
			if (jTree < 0) 
			{
				int kTree = -jTree -1;			// Need to correct for the "offset" factor

				if (kTree > locMTrees[1].size())
					cout << locTask << ", " <<  jTree << ", " << nLocHalos[1] << ", " << locMTrees[1].size() << endl;

				/* This kind of error might be due to the incorrect setting of the facRSearch variable */
				if(locMTrees[1][kTree + nLocHalos[1]].idProgenitor.size() == 0)
				{
					cout << "ERROR OnTask:" <<  locTask << ", jTree:" <<  jTree 
						<< ", nHalos:" << nLocHalos[1] << ", locTrees:" 
							<< locMTrees[1][kTree+nLocHalos[1]].progHalos.size() << endl;
	
					progHalo.Info();
					mergerTree.mainHalo.Info();
					locMTrees[1][kTree + nLocHalos[1]+1].mainHalo.Info();

				} else {
					progHalo = locMTrees[1][nLocHalos[1]+kTree].mainHalo;	
					descID = locMTrees[1][nLocHalos[1]+kTree].idProgenitor[0];	
				}
		
			} else {
				
				// Sanity check
				if (jTree > locMTrees[1].size())
					cout << "ERROR in CleanTrees(). MTree size: " << locMTrees[1].size() 
						<< ", indexj: " << jTree << endl;
						
				if (locMTrees[1][jTree].idProgenitor.size() > 0)
				{
					descID = locMTrees[1][jTree].idProgenitor[0];	
				} else {
					locMTrees[1][jTree].mainHalo.Info();
				}

				progHalo = locMTrees[1][jTree].mainHalo;
			}
#else
			descID = locMTrees[1][jTree].idProgenitor[0];	
			progHalo = locMTrees[1][jTree].mainHalo;
#endif

			/* Sanity check */
			if (descID == 0 && progHalo.nPart[1] > minPartHalo)
			{
				cout << "WARNING. Progenitor ID not assigned: " << progID << " " << descID 
					<< " | " << iTree << " " << jTree << endl;
			}

			if (mainID == descID)
			{
				mergerTree.idProgenitor.push_back(progID);
				mergerTree.indexProgenitor.push_back(jTree);
				mergerTree.progHalos.push_back(progHalo);

				for(int iT = 0; iT < nPTypes; iT++)
					mergerTree.nCommon[iT].push_back(locMTrees[0][iTree].nCommon[iT][iProg]);
			}	// mainID = descID
		}	// iProg for loop

		if (!mergerTree.isOrphan)	
			mergerTree.SortByMerit();

		if (mergerTree.idProgenitor.size() > 0)
			locCleanTrees[iStep-1].push_back(mergerTree);

		mergerTree.Clean();
	}	// kTree for loop
};



/* These two functions are used in mode = 1, when the MTrees are being read in from the .mtree files */
void AssignDescendant()
{
	uint64_t mainID = 0;
	int mainIndex = 0;

	locOrphIndex.clear();
	locOrphIndex.shrink_to_fit();
	locOrphHalos.clear();
	locOrphHalos.shrink_to_fit();

	for (int iC = 0; iC < locCleanTrees[iNumCat-1].size(); iC++)
	{
		mainID = locCleanTrees[iNumCat-1][iC].mainHalo.ID;
		mainIndex = id2Index[mainID];

		if (id2Index.find(mainID) != id2Index.end()) 
		{
			if (locCleanTrees[iNumCat-1][iC].isOrphan)
			{
				locOrphIndex.push_back(mainIndex);
				locOrphHalos.push_back(locCleanTrees[iNumCat-1][iC].mainHalo);
			}

			locCleanTrees[iNumCat-1][iC].mainHalo = locHalos[iUseCat][mainIndex];
		} else {
			cout << locTask << " does not have descendant ID: " << mainID << endl;
		}
	}

	/* Halos have been assigned, so we can clear the map */
	id2Index.clear();	
};



void AssignProgenitor()
{
	uint64_t progID = 0;
	int progIndex = 0;

	orphanHaloIndex.clear();
	orphanHaloIndex.shrink_to_fit();

	for (int iC = 0; iC < locCleanTrees[iNumCat-1].size(); iC++)
	{
		for (int iS = 0; iS < locCleanTrees[iNumCat-1][iC].progHalos.size(); iS++ )
		{
			progID = locCleanTrees[iNumCat-1][iC].progHalos[iS].ID;
	
			if(locCleanTrees[iNumCat-1][iC].isOrphan)
			{
				locCleanTrees[iNumCat-1][iC].progHalos[0] = locCleanTrees[iNumCat-1][iC].mainHalo; 
			} else {

				if (id2Index.find(progID) != id2Index.end()) 
				{
					progIndex = id2Index[progID];
#ifndef ZOOM	
					if (progIndex < 0)
						locCleanTrees[iNumCat-1][iC].progHalos[iS] = locBuffHalos[-progIndex-1];
					else
#endif
						locCleanTrees[iNumCat-1][iC].progHalos[iS] = locHalos[iUseCat][progIndex];

				} else {
					//cout << locTask << " does not have progenitor ID: " << progID << endl;
					//locCleanTrees[iNumCat-1][iC].progHalos[iS].Info();
				}
			}
		}
	}
};


/* When re-loading the merger trees initialize for each halo at z = 0 the HaloTree, which will contain all its descendants */
void InitHaloTrees()
{
	if (locTask == 0)
		cout << "Initializing halo trees..." << endl;

	locHaloTrees.resize(nLocHalos[0]);

	for (int iH = 0; iH < locHalos[0].size(); iH++) 
	{
		locHaloTrees[iH].mainHalo.resize(nSnapsUse); 
		locHaloTrees[iH].progHalo.resize(nSnapsUse);
		locHaloTrees[iH].mainHalo[0] = locHalos[0][iH]; //locCleanTrees[0][iH].mainHalo;
	}
};


/* Starting from redshift zero we build the tree backwards */
void BuildTrees()
{
	if (locTask == 0)
		cout << "Building trees..." << endl;

	map<uint64_t, int>::iterator it;

	for (int iC = 0; iC < locCleanTrees[iNumCat-1].size(); iC++)
	{
		uint64_t mainProgID = locCleanTrees[iNumCat-1][iC].progHalos[0].ID;
	
		it = id2Index.find(mainProgID);

		/* If the main progenitor is found here, otherwise add it to the list of trees to be updated */
		if (it != id2Index.end())
		{
			//locHalos[1][it->second].Info();
		} else {
				
			//cout << mainProgID << " " << it->second << " " << it->first << endl;
			// TODO Look for this halo on another task
		}
		
	
	}

};


void SyncIndex()
{
	/* Clean the map just in case */
	id2Index.clear();
	
	for (int iH = 0; iH < locHalos[iUseCat].size(); iH++)
		id2Index[locHalos[iUseCat][iH].ID] = iH;

#ifndef ZOOM
	if (iUseCat == 1)
		for (int iH = 0; iH < locBuffHalos.size(); iH++)
			id2Index[locBuffHalos[iH].ID] = -iH-1;	// Add -1 to avoid overlap with index 0
#endif
}


void DebugTrees()
{
	if (locTask == 0)
		cout << "Debugging trees for steps=" << locCleanTrees.size() << endl;
 
	for (int iC = 0; iC < locCleanTrees.size(); iC++)
	{	
		cout << "Task=(" << locTask << ") Size=(" << locCleanTrees[iC].size() << ") Step=(" << iC << ") " << endl;

		for (int iT = 0; iT < locCleanTrees[iC].size(); iT++)
				locCleanTrees[iC][iT].Info();
	}
};


