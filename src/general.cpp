#include "general.h"
#include <mpi.h>

int locTask;
int totTask;
MPI_Status status;

Halo *locHalos;
Halo *locHalosBufferSend;
Halo *locHalosBufferRecv;

size_t locHalosSize;
size_t totHalosSize;
size_t locHalosBufferSendSize;
size_t locHalosBufferRecvSize;

int nTotHalos;
int nLocHalos;

Particle *locParticles;
Particle *locParticlesBufferSend;
Particle *locParticlesBufferRecv;

size_t locParticlesSize;
size_t totParticlesSize;
size_t locParticlesBufferSendSize;
size_t locParticlesBufferRecvSize;

int nTotParticles;
int nLocParticles;

float totVmax;
float locVmax;
float bufferThickness;

