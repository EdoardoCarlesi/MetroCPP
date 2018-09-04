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
MPI_Status status;

#ifndef ZOOM
Grid GlobalGrid[2];
Grid BufferGrid;
#endif

vector<vector<Halo>> locHalos;
#ifndef ZOOM
vector<Halo> locBuffHalos;
#endif
size_t locHalosSize[2];

int iNumCat; // Halo catalog number in use, from 0 to N
int iUseCat; // Refers to 0 or 1 depending on the snapshot being used
int nTypePart = NPTYPES;
int nTotHalos[2];
int nLocHalos[2];

vector<vector<vector<vector<unsigned long long int>>>> locParts;
#ifndef ZOOM
vector<vector<vector<unsigned long long int>>> locBuffParts;
#endif
size_t locPartsSize[2];

size_t sizeHalo;
size_t sizePart;

int nPTypes;
int nLocParts[2];

/* Factors and variables */
float dMaxFactor;
float totVmax;
float locVmax;
float maxBufferThick;
float boxSize;

int nGrid;
int nChunks;
int nCat;
