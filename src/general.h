#ifndef GENERAL_H
#define GENERAL_H

#include <mpi.h>
#include <map>

#include "Halo.h"
#include "Particle.h"


class Halo;
class Particle;

// General variables and functions
unsigned int NumLines(const char *);

// MPI variables
extern int locTask;
extern int totTask;
extern MPI_Status status;

// Halo & particle related variables for each task
extern Halo *locHalos;
extern Halo *locHalosBufferSend;
extern Halo *locHalosBufferRecv;

extern size_t totHalosSize;
extern size_t locHalosSize;
extern size_t locHalosBufferSendSize;
extern size_t locHalosBufferRecvSize;

extern int nTotHalos;
extern int nLocHalos;

extern Particle **locParts;
extern Particle **locPartsBufferSend;
extern Particle **locPartsBufferRecv;

extern size_t totPartsSize;
extern size_t locPartsSize;
extern size_t locPartsBufferSendSize;
extern size_t locPartsBufferRecvSize;

extern int nTotParts;
extern int nLocParts;

// These quantities are useful to compute the buffer region
extern float totVmax;
extern float locVmax;
extern float bufferThickness;

#endif 
