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
 * global_vars.cpp: 
 * A collection of variables and structures that are being use throughout the program.
 */

#include <map>
#include <mpi.h>
#include <vector>

#include "Halo.h"
#include "Grid.h"
#include "global_vars.h"

using namespace std;

/*
 *	General variables declared here and used program-wide for synchronization
 */
int locTask;
int totTask;
MPI_Request request_recv, request_send;
MPI_Status status;

#ifndef ZOOM
Grid GlobalGrid[2];
#endif

vector<vector<MergerTree>> locMTrees;
vector<vector<MergerTree>> locCleanTrees;
vector<HaloTree> locHaloTrees;

/* We keep track of the orphan halos to synchronize them afterwards */
vector<int> orphanHaloIndex;	// Global index

vector<vector<Halo>> allHalos;
vector<vector<Halo>> locHalos;
vector<vector<vector<vector<uint64_t>>>> locParts;

#ifndef ZOOM
vector<Halo> locBuffHalos;
vector<vector<vector<uint64_t>>> locBuffParts;
#endif

vector<Halo> locOrphHalos;
vector<uint64_t> allOrphIDs;
vector<vector<vector<uint64_t>>> locOrphParts;

typedef struct Particle Particle;
vector<map <uint64_t, vector<Particle>>> locMapParts;

/* This map keeps track of the halo ids when reading from old mtree files */
vector<map<uint64_t, int>> id2Index;

map <uint64_t, int> thisMapTrees;
map <uint64_t, int> nextMapTrees;

size_t locHalosSize[2];

/* Halo catalog number in use, from 0 to N */
int iNumCat; 

/* Refers to 0 or 1 depending on the snapshot being used */
int iUseCat; 

int nTypePart = NPTYPES;

int nTotHalos[2];
int nLocHalos[2];

size_t locPartsSize[2];
size_t sizeHalo;
size_t sizePart;

int nPTypes;
int nLocParts[2];

/* Factors and variables */
float totVmax;
float locVmax;
float boxSize;

string cosmologicalModel;

int minPartCmp;
int minPartHalo;

int runMode;
int nGrid;
int facOrphanSteps;
int nTreeChunks;
int nLocChunks;
int nChunks;
int nSnapsUse;
int nSnaps;
