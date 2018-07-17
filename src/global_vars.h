#ifndef GENERAL_H
#define GENERAL_H

#include <mpi.h>
#include <map>
#include <vector>

#include "Grid.h"
#include "Halo.h"

using namespace std;


class Grid;
class Halo;

// This is the grid that keeps track of halos positions and tasks they are located on
extern Grid GlobalGrid[2];
extern Grid BufferGrid;

// MPI variables
extern int locTask;
extern int totTask;
extern MPI_Status status;

// Halo & particle related variables for each task
extern vector<vector<Halo>> locHalos;

/* Extra halos coming from the buffer nodes communicated from other tasks */
extern vector<Halo> locBuffHalos;

extern size_t locHalosSize[2];

extern int nPTypes;
extern int nTotHalos[2];
extern int nLocHalos[2];

// 0 is the last catalog (first to be read in) to the _000 one (number N)
extern int iNumCat;

// Set to 0 or 1 according to the catalogs being read in (at each step we read two at the same time)
extern int iUseCat;

/* Particles on task, also allocated by particle type within each halo */
extern vector<vector<vector<vector<unsigned long long int>>>> locParts;

/* Extra particles coming from the buffer areas located on other tasks */
extern vector<vector<vector<unsigned long long int>>> locBuffParts;

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
extern int nChunksPerFile;	// Each halo catalog / particle file is split into this number of files
#endif 
