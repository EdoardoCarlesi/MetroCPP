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
extern Grid GlobalGrid;

// MPI variables
extern int locTask;
extern int totTask;
extern MPI_Status status;

// Halo & particle related variables for each task
extern vector<Halo> locHalos;

extern size_t totHalosSize;
extern size_t locHalosSize;

extern int nPTypes;
extern int nTotHalos;
extern int nLocHalos;
extern int iLocHalos;

extern vector <vector<vector<unsigned long long int>>> locParts;	// numbers of particles are stored by particle type
extern size_t sizePart;
extern size_t sizeHalo;

/* Particle and halo sizes are communicated to check for integrity */
extern size_t totPartsSize;
extern size_t locPartsSize;

extern int nTotParts;
extern int nLocParts;

// These quantities are useful to compute the buffer region
extern float locXmin[3];
extern float locXmax[3];

extern float totVmax;
extern float locVmax;
extern float bufferThickness;
extern int nChunksPerFile;	// Each halo catalog / particle file is split into this number of files

#endif 
