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

// General variables and functions
unsigned int NumLines(const char *);
void InitLocVariables(void);
float VectorModule(float *);
void CleanMemory(void);

// This is the grid that keeps track of halos positions and tasks they are located on
extern Grid GlobalGrid[2];

// MPI variables
extern int locTask;
extern int totTask;
extern MPI_Status status;

// Halo & particle related variables for each task
extern vector<vector<Halo>> locHalos;

extern size_t locHalosSize[2];

extern int nPTypes;
extern int nTotHalos[2];
extern int nLocHalos[2];

// 0 is the last catalog (first to be read in) to the _000 one (number N)
extern int iNumCat;
// Set to 0 or 1 according to the catalogs being read in
extern int iUseCat;

// numbers of particles are stored by particle type
extern vector<vector<vector<vector<unsigned long long int>>>> locParts;

extern size_t sizePart;
extern size_t sizeHalo;

/* Particle and halo sizes are communicated to check for integrity */
extern size_t locPartsSize[2];
extern int nLocParts[2];

extern float totVmax;
extern float locVmax;
extern float bufferThickness;
extern int nChunksPerFile;	// Each halo catalog / particle file is split into this number of files
#endif 
