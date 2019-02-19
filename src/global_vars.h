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


#ifndef GLOBAL_VARS_H
#define GLOBAL_VARS_H

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

/* HaloTree contains a mainHalo at z=0 and its full history */
extern vector<HaloTree> locHaloTrees;

/* Here we store the pairwise connections between halos in catalog 0 and 1, before being sorted */
extern vector<vector<MergerTree>> locMTrees;

/* This variables stores the number of (clean) connections between halos in catalog 0 and catalog 1, two steps at the time */
extern vector<vector<MergerTree>> locCleanTrees;

/* Helps connecting halos when rebuilidng the trees from input files */
extern vector<map<uint64_t, int>> id2Index;

/* Particles on task, also allocated by particle type within each halo */
extern vector<vector<vector<vector<uint64_t>>>> locParts;

#ifndef ZOOM
/* Extra halos coming from the buffer nodes communicated from other tasks */
extern vector<Halo> locBuffHalos;
extern vector<vector<vector<uint64_t>>> locBuffParts;
#endif

/* These vectors keep track of orphan halos for which no progenitor could be found (so far) */
extern vector<uint64_t> allOrphIDs;
extern vector<Halo> locOrphHalos;
extern vector<vector<vector<uint64_t>>> locOrphParts;

struct Particle {
	int type;
	uint64_t haloID;
};

extern vector<map<uint64_t, vector<Particle>>> locMapParts;

extern map <uint64_t, int> thisMapTrees;
extern map <uint64_t, int> nextMapTrees;

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

extern float boxSize;
extern float totVmax;
extern float locVmax;

extern string cosmologicalModel;

/* These int values are being read from the configuration file */
extern int runMode;

// Track orphan halos to a maximum of these steps
extern int facOrphanSteps;	

extern int nSnapsUse;
extern int nSnaps;

// Number of MPI tasks used when writing the tree
extern int nTreeChunks;	

// Each halo catalog / particle file is split into this number of files
extern int nChunks;	

// Each tast has a local number of chunks to read (it should be equal for all tasks for better load balancing, but in general it can vary)
extern int nLocChunks;  
#endif 
