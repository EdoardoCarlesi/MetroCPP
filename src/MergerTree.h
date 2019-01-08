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

	vector<Halo> progHalos;				// progenitor Halos of the main tree
	Halo mainHalo;					// Main halo of the MTree, to be stored in the cleantree only

	bool isOrphan;					// If no progenitor is found, the halo is orphan and a token placeholder halo is 
							// created with the same particle content to keep tracking it at subsequent steps
	vector<unsigned long long int> idProgenitor;	// IDs of progenitors
	vector<int> indexProgenitor;			// local array index of progenitors
	vector<vector<int>> nCommon;			// Particles in common are separated per particle type

#ifdef CMP_MAP
	map<unsigned long long int, vector<int>> indexCommon; 	
	void AssignMap(void);
#endif

	void SortByMerit(void);				// Once possible progenitors have been found, compare
	void Clean(void);
	void Info(void);
};


/* This class stores the main halo and its progenitors at each step */
class HaloTree {
	
public:
	HaloTree();
	~HaloTree();

	int nStep;
	
	vector<MergerTree> mTree; 
	vector<Halo> mainHalo;				// This traces the main progenitor

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

void AssignDescendant(void);
void AssignProgenitor(void);

// Pairwise comparison of halos
void FindProgenitors(int, int);

// Decide whether to compare two halos
bool CompareHalos(int, int, int, int);

// Given two (sorted) vectors, compare their content and return the number of common elements
vector<int> CommonParticles(vector<vector<unsigned long long int>>, vector<vector<unsigned long long int>>);

#endif
