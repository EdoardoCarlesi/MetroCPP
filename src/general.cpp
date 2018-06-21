#include "general.h"
#include <mpi.h>

#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

using namespace std;


int locTask;
int totTask;
MPI_Status status;

vector <Halo> locHalos;
vector <Halo> locHalosBufferSend;
vector <Halo> locHalosBufferRecv;

size_t locHalosSize;
size_t totHalosSize;
size_t locHalosBufferSendSize;
size_t locHalosBufferRecvSize;

int nTotHalos;
int nLocHalos;

vector <vector <Particle>> locParts;
void *locPartsBufferSend;
void *locPartsBufferRecv;

size_t sizeHalo;
size_t sizePart;
size_t locPartsSize;
size_t totPartsSize;
size_t locPartsBufferSendSize;
size_t locPartsBufferRecvSize;

int nTotParts;
int nLocParts;

// Apex coordinates of the sub-box containing the halo on each task
float locXmin[3];
float locXmax[3];

float totVmax;
float locVmax;
float bufferThickness;
int nChunksPerFile;


/* 
 * 	General functions
 */

void InitLocVariables(void)
{
        locVmax = 0.0; 
	totVmax = 0.0;
        locHalosSize = 0;
        totHalosSize = 0;

	for (int iX = 0; iX < 3; iX++)
	{
		locXmin[iX] = 10000000.0;
		locXmax[iX] = 0.0;
	}

	sizeHalo = sizeof(Halo);
        sizePart = sizeof(Particle);
	
	nTotParts = 0; nLocParts = 0;
};


unsigned int NumLines(const char * fileName)
{
	unsigned int nLines = 0;
	const char * charHead = "#";

	string lineIn;

	ifstream fileIn(fileName);

	if (fileIn)
	{
		while (getline(fileIn, lineIn))
		{
			const char * lineTmp = lineIn.c_str();

			if (lineTmp[0] != charHead[0])
				nLines++;
		}

		return nLines;
	}
	else
	{ 
		cout << "File " <<  fileName << " not found on task: " << locTask << endl;
		return -1;
	}
};


float VectorModule(float *V)
{
	return sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
};

