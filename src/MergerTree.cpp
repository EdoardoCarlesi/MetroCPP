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

#include <algorithm>
#include <string>
#include <vector>
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


/* Append a list of progenitors from another tree to this one */
void MergerTree::Append(MergerTree mTree)
{
	if (mainHalo.ID != mTree.mainHalo.ID)
		cout << "WARNING. Appending mtree of halo " << mTree.mainHalo.ID << " to a different main branch: " << mainHalo.ID << endl;

	for (int iM = 0; iM < mTree.progHalo.size(); iM ++)
	{
		map<uint64_t, vector<int>>::iterator iter;
		iter = indexCommon.find(mTree.progHalo[iM].ID);

		/* If the progenitor is not inside this tree then append it */
		if (iter == indexCommon.end())
		{
			idProgenitor.push_back(mTree.progHalo[iM].ID);
			progHalo.push_back(mTree.progHalo[iM]);

			for (int iC = 0; iC < nPTypes; iC ++)
				nCommon[iC].push_back(mTree.nCommon[iC][iM]);
		}
	}

};


void MergerTree::Info()
{
	cout << "Task=" << locTask << " " << mainHalo.ID << " " << idProgenitor.size() << " nPart: " << mainHalo.nPart[1] << endl;

	/*
	if (isOrphan)
		for (int iP = 0; iP < idProgenitor.size(); iP++)
			cout << " ID= " << idProgenitor[iP] << endl;
	else
	*/
	for (int iP = 0; iP < idProgenitor.size(); iP++)
		cout << iP << " ID= " << idProgenitor[iP] << " nPartComm=" << nCommon[1][iP] <<  " nPart: " << progHalo[iP].nPart[1] << endl;
};



void MergerTree::Clean()
{
	for (int iC = 0; iC < nCommon.size(); iC++)
	{
		nCommon[iC].clear();
		nCommon[iC].shrink_to_fit();
	}

	nCommon.clear();
	nCommon.shrink_to_fit();

	idProgenitor.clear();
	idProgenitor.shrink_to_fit();
	progHalo.clear();
	progHalo.shrink_to_fit();
};



MergerTree::MergerTree()
{
	nCommon.resize(nPTypes);
};

MergerTree::~MergerTree()
{
};


void MergerTree::AssignMap()
{
	int nProgs = 0, nCommTot = 0, iP = 0;
	
	/* Each merger tree stores the halo ids in a map that stores the number of particles shared with each progenitor */
	for (auto const& thisMap : indexCommon) 
	{
		uint64_t thisHaloID = thisMap.first;
		vector<int> thisCommon = thisMap.second;

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
		//cout << "Halo " << mainHalo.ID << " has no progenitors. " << endl;
	} else { 
		/* We allocate these vectors here, the indexes and the halos will be copied inside later on */
		isOrphan = false;
		progHalo.resize(nProgs);
	}
};


/* This sorting algorithm might be very inefficient but it's straightforward to implement, 
 * plus we will rarely deal with halos with more than 10^3 progenitors */
void MergerTree::SortByMerit()
{
	vector<uint64_t> tmpIdx;
	vector<vector<int>> tmpNCommon;
	vector<float> allMerit;
	vector<int> idx;
	vector<Halo> tmpProgHalo;
	float merit = 0.0;

	for (int iM = 0; iM < progHalo.size(); iM++)
	{
		int nComm = 0, nPH0 = 0, nPH1 = 0;
		double ratioM = 0;

		for (int iC = 0; iC < nPTypes; iC++)
			nComm += nCommon[iC][iM];
 
		for (int iP = 0; iP < nPTypes; iP++)
		{	
			nPH0 += mainHalo.nPart[iP];
			nPH1 += progHalo[iM].nPart[iP];
		}

		ratioM = (float) nPH0 / (float) nPH1;
	
		if (ratioM < 1.0) ratioM = 1.0 / ratioM;

		merit = nComm  / (ratioM*1.01 - 1.0);
		merit *= (1.0 + 0.00001 * iM);	// We change the merit slightly, 
						// to avoid confusion when two halos have the same number of particles 
						// and the same number of particles shared with the host halo

		allMerit.push_back(merit);
	}
	
	idx = SortIndexes(allMerit);
	tmpProgHalo.resize(idx.size());
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
		tmpProgHalo[oldIdx] = progHalo[newIdx];
	}	

	for (int iM = 0; iM < idProgenitor.size(); iM++)
	{
		idProgenitor[iM] = tmpIdx[iM];
		progHalo[iM] = tmpProgHalo[iM];

		for (int iT = 0; iT < nPTypes; iT++)
			nCommon[iT][iM] = tmpNCommon[iT][iM];
	}

	tmpIdx.clear();
	tmpIdx.shrink_to_fit();
	tmpProgHalo.clear();
	tmpProgHalo.shrink_to_fit();
};

	       /****************************************************************************
		* These are general functions which do not belong to the class Merger Tree *
		****************************************************************************/

#ifndef ZOOM		

/* This is a very fast way of comparing particles content of halos across snapshots, that relies on maps 
 * and scales linearly with the number of particles on each task */

void FindProgenitors(int iOne, int iTwo)
{
	int nLoopHalos[2], iOldOrphans = 0, iFixOrphans = 0, nLocOrphans = 0, nLocTrees = 0, nLocUntrack = 0; 

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

			if (thisHaloIndex >= nLocHalos[iTwo])
				locMTrees[iOne][iM].progHalo[iP] = locBuffHalos[thisHaloIndex-nLocHalos[iTwo]];
			else
				locMTrees[iOne][iM].progHalo[iP] = locHalos[iTwo][thisHaloIndex];
		}

		/* Sort only with 2 progenitors at least */
		if (locMTrees[iOne][iM].progHalo.size() > 1)
			locMTrees[iOne][iM].SortByMerit();
	}
};		/* End of the find progenitor function in full box mode */

#else		 /* The Find Progenitors function in ZOOM MODE is different than the standard one */

void FindProgenitors(int iOne, int iTwo)
{
	int iOldOrphans = 0, iFixOrphans = 0, nLocOrphans = 0, nLocTrees = 0, nLocUntrack = 0; 

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
		int nProgs = locMTrees[iOne][iM].idProgenitor.size();
		int thisHaloIndex = 0;
		uint64_t thisHaloID;

		for (int iP = 0; iP < nProgs; iP++)
		{
			thisHaloID = locMTrees[iOne][iM].idProgenitor[iP];
			thisHaloIndex = nextMapTrees[thisHaloID];			
			locMTrees[iOne][iM].progHalo[iP] = locHalos[iTwo][thisHaloIndex];
		}

		/* Sort only with 2 progenitors at least */
		if (locMTrees[iOne][iM].progHalo.size() > 1)
			locMTrees[iOne][iM].SortByMerit();
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
	int thisIndex = 0, nLocUntrack = 0, nLocOrphans = 0; 

	if (locTask == 0)
		cout << "Cleaning Merger Tree connections for " << locMTrees[0].size() << " halos." << endl;

#ifdef VERBOSE
	if (locTask == 0)
		cout << "nHalos " << locHalos[0].size() << ", nTrees: " << locMTrees[0].size() << endl; 
#endif

	for (int iTree = 0; iTree < locMTrees[0].size(); iTree++)
	{
		uint64_t mainID = locMTrees[0][iTree].mainHalo.ID;
		int nProgSize = locMTrees[0][iTree].idProgenitor.size();

		MergerTree mergerTree;
		mergerTree.mainHalo = locMTrees[0][iTree].mainHalo;

		/* At each step we only record the connections between halos in catalog 0 and catalog 1, without attempting at a
		 * reconstruction of the full merger history. This will be done later. */
		for (int iProg = 0; iProg < nProgSize; iProg++)
		{
			Halo progHalo;	
			uint64_t progID = locMTrees[0][iTree].idProgenitor[iProg];
			uint64_t descID;
	
			/* 1-trees are on thisMap while 0-trees are on nextMap */
			int jTree = thisMapTrees[progID]; 

			if (jTree > locMTrees[1].size() || jTree < 0) 
			{
				cout << " ON TASK " << locTask << " jTree is outside the limits: " << jTree << endl; 
			} else {
				descID = locMTrees[1][jTree].progHalo[0].ID;	
				progHalo = locMTrees[1][jTree].mainHalo;
			}

			/* Sanity check */
			if (descID == 0 && progHalo.nPart[1] > minPartHalo)
			{
				//locMTrees[1][jTree].Info(); 
				cout << "OnTask= " << locTask << ": WARNING, progen. ID: " << progID << " has no descID: " << descID 
					<< " | " << iTree << " " << jTree << endl;
			}

			/* Check that the main halo is the likeliest descendant of its progenitor */
			if (mainID == descID)
			{
				//cout << "main: " << mainID << " desc: " << descID << " , j:" << jTree << " , i:" << iTree << endl; 
				mergerTree.idProgenitor.push_back(progID);
				mergerTree.progHalo.push_back(progHalo);

				for(int iT = 0; iT < nPTypes; iT++)
					mergerTree.nCommon[iT].push_back(locMTrees[0][iTree].nCommon[iT][iProg]);
			}	// mainID = descID
		}	// loop on progenitors

		/* In some cases, a subhalo disappears and its main progenitor turns out to be the host halo, 
		 * In this case, the subhalo is not recorded among the orphan halos, since it does have a connection
		 * and shared particles in the forward loop. Here we check again that this subhalo is not the main 
		 * descendent of a progenitor host, and record it among the orphan halos to be tracked */
		if (mergerTree.idProgenitor.size() == 0 && mergerTree.mainHalo.nAllPart() > minPartHalo) 
		{
			Halo thisHalo = mergerTree.mainHalo;
			thisHalo.nOrphanSteps++;
			thisHalo.isToken = true;

			int locMaxOrphanSteps = 1 + int (thisHalo.nAllPart() / facOrphanSteps);

			/* Upper limit on the total number of steps an halo can be tracked */
			if (locMaxOrphanSteps > maxOrphanSteps)		
				locMaxOrphanSteps = maxOrphanSteps;

			/* Check if it's worth to continue tracking this orphan halo */
			if (thisHalo.nOrphanSteps <= locMaxOrphanSteps)
			{
#ifdef GATHER_TREES
				allOrphIDs.push_back(thisHalo.ID);
#else

#ifdef ZOOM
#endif
				/* Update the container of local orphan halos */
				locOrphHalos.push_back(thisHalo);

				/* Update the particle content */
				locOrphParts.push_back(locParts[0][iTree]);
				locOrphParts[nLocOrphans].resize(nPTypes);

				/* Copy particle blocks divided by particle type */
				for (int iP = 0; iP < nPTypes; iP++)
					copy(locParts[0][iTree][iP].begin(), locParts[0][iTree][iP].end(), 
						back_inserter(locOrphParts[nLocOrphans][iP]));
#ifdef ZOOM
#endif
#endif
				mergerTree.isOrphan = true;
				mergerTree.idProgenitor.push_back(thisHalo.ID);
				mergerTree.progHalo.push_back(thisHalo);

				for(int iT = 0; iT < nPTypes; iT++)
					mergerTree.nCommon[iT].push_back(locMTrees[0][iTree].mainHalo.nPart[iT]);

				nLocOrphans++;

			} else {	// We stop following this orphan halo, too small and disconnected for too many steps
				nLocUntrack++;
			}	
 
		} else {	// if the tree has no progenitors  
			mergerTree.isOrphan = false;
		}	

		if (mergerTree.idProgenitor.size() > 1)	
			mergerTree.SortByMerit();

		if (mergerTree.idProgenitor.size() > 0) 
			locCleanTrees[iStep-1].push_back(mergerTree);

		mergerTree.Clean();
	}	// iTree for loop

	/* Final statistics - sanity check */
	int nTotOrphans = 0, nTotUntrack = 0;

#ifdef GATHER_TREES
	nLocOrphans = allOrphIDs.size();
	nTotOrphans = nLocOrphans;
	nTotUntrack = nLocUntrack;
#else
	nLocOrphans = locOrphHalos.size();
	MPI_Reduce(&nLocOrphans,   &nTotOrphans, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&nLocUntrack,   &nTotUntrack, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

#endif

	if (locTask == 0)
		cout << "The total number of orphan halos after cleaning the connections is: " << nTotOrphans
			<< " while " << nTotUntrack << " will be untracked. "  << endl;
};



/* These two functions are used in mode = 1, when the MTrees are being read in from the .mtree files */
void AssignDescendant()
{
	uint64_t mainID = 0;
	int mainIndex = 0;

	locOrphHalos.clear();
	locOrphHalos.shrink_to_fit();

	for (int iC = 0; iC < locCleanTrees[iNumCat-1].size(); iC++)
	{
		mainID = locCleanTrees[iNumCat-1][iC].mainHalo.ID;
		mainIndex = id2Index[0][mainID];

		if (id2Index[0].find(mainID) != id2Index[0].end()) 
		{
			if (locCleanTrees[iNumCat-1][iC].isOrphan)
			{
				locOrphHalos.push_back(locCleanTrees[iNumCat-1][iC].mainHalo);
			}

			locCleanTrees[iNumCat-1][iC].mainHalo = locHalos[iUseCat][mainIndex];
		} else {
			cout << locTask << " does not have descendant ID: " << mainID << endl;
		}
	}
};



void AssignProgenitor()
{
	uint64_t progID = 0;
	int progIndex = 0;

	for (int iC = 0; iC < locCleanTrees[iNumCat-1].size(); iC++)
	{
		for (int iS = 0; iS < locCleanTrees[iNumCat-1][iC].progHalo.size(); iS++ )
		{
			progID = locCleanTrees[iNumCat-1][iC].progHalo[iS].ID;
	
			if(locCleanTrees[iNumCat-1][iC].isOrphan)
			{
				locCleanTrees[iNumCat-1][iC].progHalo[0] = locCleanTrees[iNumCat-1][iC].mainHalo; 
			} else {

				if (id2Index[1].find(progID) != id2Index[1].end()) 
				{
					progIndex = id2Index[1][progID];
#ifndef ZOOM	
					if (progIndex < 0)
						locCleanTrees[iNumCat-1][iC].progHalo[iS] = locBuffHalos[-progIndex-1];
					else
#endif
						locCleanTrees[iNumCat-1][iC].progHalo[iS] = locHalos[iUseCat][progIndex];

				} else {
					//cout << locTask << " does not have progenitor ID: " << progID << endl;
					//locCleanTrees[iNumCat-1][iC].progHalo[iS].Info();
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

	id2Index.resize(2);
	locHaloTrees.resize(nLocHalos[0]);

	for (int iH = 0; iH < locCleanTrees[0].size(); iH++) 
	{
		locHaloTrees[iH].mainHalo.resize(nSnapsUse); 
		locHaloTrees[iH].progHalo.resize(nSnapsUse);
		locHaloTrees[iH].mainHalo[0] = locCleanTrees[0][iH].mainHalo;

		//if (locCleanTrees[0][iH].mainHalo.ID != locHalos[0][iH].ID)
		//	locHalos[0][iH].Info();
			//cout << iH << " " << locHalos[0][iH].ID << endl;
	}
};


/* Starting from redshift zero we build the tree backwards */
void BuildTrees()
{
	if (locTask == 0)
		cout << "Building trees..." << endl;

	map<uint64_t, int>::iterator it;

	int iFound = 0, iNotFound = 0, iOrph = 0;
	int nHaloTrees = locCleanTrees[iNumCat-1].size();

	for (int iC = 0; iC < nHaloTrees; iC++)
	{
		uint64_t mainProgID = locCleanTrees[iNumCat-1][iC].progHalo[0].ID;
	
		it = id2Index[1].find(mainProgID);

		/* If the main progenitor is found here, otherwise add it to the list of trees to be updated */
		if (it != id2Index[1].end())
		{
			//if (it->second > nLocHalos[1])
			//	cout << "buffer halo: " << it-> second << endl;

			iFound++;	
		} else {

			if (locCleanTrees[iNumCat-1][iC].isOrphan)
			{
				iOrph++;
				// It's an orphan halo, just replicate it at the next step
			} else {
				// If this is not an orphan halo, we may need to look for it on another task
				iNotFound++;
			}
				/*
					cout << iC << ", Orphan halo " << locCleanTrees[iNumCat-1][iC].mainHalo.ID << endl;
					if (locCleanTrees[iNumCat-1][iC].progHalos[0].isToken);
					cout << iC << ", Token halo " << locCleanTrees[iNumCat-1][iC].progHalos[0].isToken << endl;
				*/

			//locCleanTrees[iNumCat-1][iC].mainHalo.Info();
			//locCleanTrees[iNumCat-1][iC].progHalos[0].Info();
			//cout << mainProgID << " " << it->second << " " << it->first << endl;
		}
	}
	
	if (locTask == 0)
		cout << "OnTask = " << locTask << " nHalos: " << nHaloTrees << " found: " << iFound 
			<< " not found: " << iNotFound << " orphans: " << iOrph << endl; 

};


void SyncIndex()
{
	/* Clean the map just in case */
	id2Index[iUseCat].clear();
	
	for (int iH = 0; iH < locHalos[iUseCat].size(); iH++)
		id2Index[iUseCat][locHalos[iUseCat][iH].ID] = iH;

#ifndef ZOOM
	if (iUseCat == 1)
	{
		//cout << "SyncIndex for buffer halos: " << locBuffHalos.size() << endl;
 
		for (int iH = 0; iH < locBuffHalos.size(); iH++)
			id2Index[1][locBuffHalos[iH].ID] = nLocHalos[1] + iH;	
	}
#endif
}


void FreeMergerTrees(int iNumCat)
{
	if (locTask == 0)
		cout << "Freeing MergerTrees... " << endl;

        for (auto thisMTree : locMTrees[0])
                thisMTree.Clean();

        locMTrees[0].clear();
        locMTrees[0].shrink_to_fit();

        for (auto thisMTree : locMTrees[1])
                thisMTree.Clean();

        locMTrees[1].clear();
        locMTrees[1].shrink_to_fit();

	for (auto thisMTree : locCleanTrees[iNumCat-1])
		thisMTree.Clean();

	locCleanTrees[iNumCat-1].clear();
	locCleanTrees[iNumCat-1].shrink_to_fit();

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


