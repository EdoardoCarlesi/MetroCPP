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


#ifndef GENERAL_H
#define GENERAL_H

#include <mpi.h>
#include <map>
#include <vector>

#include "Grid.h"
#include "Halo.h"
#include "MergerTree.h"

using namespace std;

class MergerTree;
class Halo;

#ifndef ZOOM
class Grid;

// This is the grid that keeps track of halos positions and tasks they are located on
extern Grid GlobalGrid[2];
#endif

/* MPI variables */
extern int locTask;
extern int totTask;
extern MPI_Status status;
extern MPI_Request request_recv, request_send;

/* Halo & particle related variables for each task */
extern vector<vector<Halo>> locHalos;

/* Here we store all halo catalogs at all redshifts */
extern vector<vector<Halo>> allHalos;

/* Here we store the pairwise connections between halos in catalog 0 and 1, before being sorted */
extern vector<vector<MergerTree>> locMTrees;

/* This variables stores the number of (clean) connections between halos in catalog 0 and catalog 1, two steps at the time */
extern vector<vector<MergerTree>> locCleanTrees;

/* Helps connecting halos when rebuilidng the trees from input files */
extern map <unsigned long long int, int> id2Index;

/* Particles on task, also allocated by particle type within each halo */
extern vector<vector<vector<vector<unsigned long long int>>>> locParts;

#ifndef ZOOM
/* Extra halos coming from the buffer nodes communicated from other tasks */
extern vector<Halo> locBuffHalos;
extern vector<vector<vector<unsigned long long int>>> locBuffParts;
#else

/* In full box mode each task keeps track of its orphans locally */
extern vector<Halo> locOrphHalos;
extern vector<int> locOrphIndex;

/* Extra particles coming from the buffer areas located on other tasks */
extern vector<vector<vector<unsigned long long int>>> locOrphParts;
#endif

struct Particle {
	int type;
	unsigned long long int haloID;
};

extern vector<map<unsigned long long int, vector<Particle>>> locMapParts;

extern map <unsigned long long int, int> thisMapTrees;
extern map <unsigned long long int, int> nextMapTrees;

extern size_t locHalosSize[2];

extern int minPartCmp;
extern int minPartHalo;

// We keep track of the orphan halos to synchronize them afterwards
extern vector<int> orphanHaloIndex;

extern int nPTypes;
extern int nTotHalos[2];
extern int nLocHalos[2];

// 0 is the last catalog (first to be read in) to the _000 one (number N)
extern int iNumCat;

// Set to 0 or 1 according to the catalogs being read in (at each step we read two at the same time)
extern int iUseCat;

extern size_t sizePart;
extern size_t sizeHalo;

/* Particle and halo sizes are communicated to check for integrity */
extern size_t locPartsSize[2];
extern int nLocParts[2];
extern int nGrid;

extern float dMaxFactor;
extern float boxSize;
extern float totVmax;
extern float locVmax;
extern float maxBufferThick;

extern string cosmologicalModel;

/* These int values are being read from the configuration file */
extern int runMode;
extern int facOrphanSteps;	// Track orphan halos to a maximum of these steps
extern int nSnapsUse;
extern int nSnaps;
extern int nTreeChunks;	// Number of MPI tasks used when writing the tree
extern int nChunks;	// Each halo catalog / particle file is split into this number of files
extern int nLocChunks;  // Each tast has a local number of chunks to read (it should be equal for all tasks for better load balancing, but in general it can vary)
#endif 
