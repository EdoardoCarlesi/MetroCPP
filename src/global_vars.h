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
extern Grid BufferGrid;
#endif

/* MPI variables */
extern int locTask;
extern int totTask;
extern MPI_Status status;

/* Halo & particle related variables for each task */
extern vector<vector<Halo>> locHalos;

/* Here we store all halo catalogs at all redshifts */
extern vector<vector<Halo>> allHalos;

/* Here we store the pairwise connections between halos in catalog 0 and 1, before being sorted */
extern vector<vector<MergerTree>> locMTrees;

#ifdef ZOOM
/* This variable keeps track of the halos that should be used on each task for the backward comparison, 
 * to avoid looping on all the halos */
#endif
extern vector<int> locTreeIndex;	

/* This variables stores the number of (clean) connections between halos in catalog 0 and catalog 1, two steps at the time */
extern vector<vector<MergerTree>> locCleanTrees;

/* Maps that contain halo ids and a link to the locHalo index */
#ifndef ZOOM
extern map <unsigned long long int, int> locId2Index;
#endif
extern map <unsigned long long int, int> id2Index;

#ifndef ZOOM
/* Extra halos coming from the buffer nodes communicated from other tasks */
extern vector<Halo> locBuffHalos;
#endif

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

/* Particles on task, also allocated by particle type within each halo */
extern vector<vector<vector<vector<unsigned long long int>>>> locParts;

#ifndef ZOOM
/* Extra particles coming from the buffer areas located on other tasks */
extern vector<vector<vector<unsigned long long int>>> locBuffParts;
#endif

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
#endif 
