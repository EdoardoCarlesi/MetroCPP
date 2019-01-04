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

extern map <unsigned long long int, int> id2Index;

/* Particles on task, also allocated by particle type within each halo */
extern vector<vector<vector<vector<unsigned long long int>>>> locParts;

#ifndef ZOOM
/* Map that contain halo ids and a link to the locHalo index */
extern map <unsigned long long int, int> locId2Index;

/* Extra halos coming from the buffer nodes communicated from other tasks */
extern vector<Halo> locBuffHalos;

/* In full box mode each task keeps track of its orphans locally */
extern vector<Halo> locOrphHalos;
extern vector<int> locOrphIndex;

/* Extra particles coming from the buffer areas located on other tasks */
extern vector<vector<vector<unsigned long long int>>> locBuffParts;

extern vector<vector<vector<unsigned long long int>>> locOrphParts;
#else

/* This variable keeps track of the halos that should be used on each task for the backward comparison, 
 * to avoid looping on all the halos */
extern vector<int> locTreeIndex;	
#endif

#ifdef CMP_MAP	// The key is the particle ID, the result is a vector of Halo IDs it belongs to
struct Particle {
	int type;
	unsigned long long int haloID;
};

extern vector<map<unsigned long long int, vector<Particle>>> locMapParts;
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
