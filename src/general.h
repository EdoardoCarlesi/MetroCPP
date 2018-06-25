#ifndef GENERAL_H
#define GENERAL_H

#include <mpi.h>
#include <map>
#include <vector>

#include "Grid.h"
#include "Halo.h"
#include "Particle.h"

using namespace std;

class Grid;
class Halo;
class Particle;

// General variables and functions
void InitLocVariables(void);
unsigned int NumLines(const char *);
float VectorModule(float *);

extern Grid GlobalGrid;

// MPI variables
extern int locTask;
extern int totTask;
extern MPI_Status status;

// Halo & particle related variables for each task
extern vector<Halo> locHalos;
extern vector<Halo> locHalosBufferSend;
extern vector<Halo> locHalosBufferRecv;

extern size_t totHalosSize;
extern size_t locHalosSize;
extern size_t locHalosBufferSendSize;
extern size_t locHalosBufferRecvSize;

extern int nTotHalos;
extern int nLocHalos;

extern vector<vector <Particle>> locParts;
extern vector<vector <Particle>> locPartsBufferSend;
extern vector<vector <Particle>> locPartsBufferRecv;

extern size_t sizePart;
extern size_t sizeHalo;

/* Particle and halo sizes are communicated to check for integrity */
extern size_t totPartsSize;
extern size_t locPartsSize;
extern size_t locPartsBufferSendSize;
extern size_t locPartsBufferRecvSize;

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
