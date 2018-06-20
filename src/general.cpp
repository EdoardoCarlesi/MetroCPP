#include "general.h"
#include <mpi.h>

#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;


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

Particle **locParts;
Particle **locPartsBufferSend;
Particle **locPartsBufferRecv;

size_t locPartsSize;
size_t totPartsSize;
size_t locPartsBufferSendSize;
size_t locPartsBufferRecvSize;

int nTotParts;
int nLocParts;

float totVmax;
float locVmax;
float bufferThickness;


/* 
 * 	General functions
 */

unsigned int NumLines(const char * fileName)
{
	unsigned int nLines = 0;
	string lineIn;

	ifstream fileIn(fileName);

	if (fileIn)
	{
		while (getline(fileIn, lineIn))
			nLines++;

		return nLines;
	}
	else
	{ 
		cout << "File " <<  fileName << " not found on task: " << locTask << endl;
		return -1;
	}
};

