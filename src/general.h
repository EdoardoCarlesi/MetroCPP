#ifndef GENERAL_H
#define GENERAL_H

#include <mpi.h>
#include <map>

#include "Halo.h"
#include "Particle.h"

//typedef unsigned long long int 	ullint;
//typedef std::map<ullint, unsigned int> idIndex;

class Halo;
class Particle;

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

extern Particle *locParticles;
extern Particle *locParticlesBufferSend;
extern Particle *locParticlesBufferRecv;

extern size_t totParticlesSize;
extern size_t locParticlesSize;
extern size_t locParticlesBufferSendSize;
extern size_t locParticlesBufferRecvSize;

extern int nTotParticles;
extern int nLocParticles;

// These quantities are useful to compute the buffer region
extern float totVmax;
extern float locVmax;
extern float bufferThickness;

#endif 
