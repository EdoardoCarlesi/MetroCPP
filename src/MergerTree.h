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


#ifndef MERGERTREE_H
#define MERGERTREE_H

#include <map>
#include <string>

#include "Halo.h"

using namespace std;


/* Due to the reliance of this class on vector template, we cannot directly MPI_Sendrecv the MergerTrees */
class MergerTree {

public:
	MergerTree();
	~MergerTree();

	Halo mainHalo;					// Main halo of the MTree, to be stored in the cleantree only
	vector<Halo> progHalo;				// progenitor Halos of the main tree

	bool isOrphan;					// If no progenitor is found, the halo is orphan and a token placeholder halo is 
							// created with the same particle content to keep tracking it at subsequent steps
	vector<uint64_t> idProgenitor;			// IDs of progenitors --> this is needed to track halos from maps and then load the progenitors
	vector<vector<int>> nCommon;			// Particles in common are separated per particle type

	map<uint64_t, vector<int>> indexCommon; 	

	void Append(MergerTree);
	void SortByMerit(void);				// Once possible progenitors have been found, compare
	void AssignMap(void);				
	void Clean(void);
	void Info(void);
};


/* This class stores the main halo and its progenitors at each step */
class HaloTree {
	
public:
	HaloTree();
	~HaloTree();

	int nStep;
	
	/* main descendant halo at z=0 */
	vector<Halo> mainHalo;				

	/* Vector of n steps, each step containing all progenitor halos */
	vector<vector<Halo>> progHalo;


	void SmoothTree(void);				// Smooths over fly-bys 
	void FixTree(void);				// Looks for missing subhalos and fixes with token halos at the missing positions
	void Clean(void);

	void WriteMergerTree(void);			// Prints all the informations 
	void WriteTrajectory(void);
	void WriteMAH(void);
	void WriteIDs(void);				// Only print the ID of each halo and its main progenitor
};


/*
 * Functions used to build the merger trees 
 */

void InitTrees(int);
void CleanTrees(int);
void DebugTrees(void);

/* These functions are used in mode 1, when reading in a raw set of merger tree files */
void AssignDescendant(void);
void AssignProgenitor(void);
void InitHaloTrees(void);
void SyncIndex(void);
void BuildTrees(void);

// Pairwise comparison of halos
void FindProgenitors(int, int);

#ifdef ZOOM
// Decide whether to compare two halos
bool CompareHalos(int, int, int, int);

// Given two (sorted) vectors, compare their content and return the number of common elements
vector<int> CommonParticles(vector<vector<uint64_t>>, vector<vector<uint64_t>>);
#endif

#endif
