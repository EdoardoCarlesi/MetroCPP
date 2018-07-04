#include <mpi.h>

#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include "Grid.h"
#include "general.h"

using namespace std;


int locTask;
int totTask;
MPI_Status status;

Grid GlobalGrid[2];

vector<vector<Halo>> locHalos;
vector<Halo> locHalosBuffer;
size_t locHalosSize[2];

int iNumCat; // Halo catalog number in use, from 0 to N
int iUseCat; // Refers to 0 or 1 depending on the snapshot being used
int nTypePart;
int nTotHalos[2];
int nLocHalos[2];

vector<vector<vector<vector<unsigned long long int>>>> locParts;
vector<vector<vector<unsigned long long int>>> locPartsBuffer;
size_t locPartsSize[2];

size_t sizeHalo;
size_t sizePart;

int nPTypes;
int nLocParts[2];

float totVmax;
float locVmax;
float maxBufferThick;
float boxSize;

int nGrid;
int nChunksPerFile;


/* 
 * 	General functions
 */
void InitLocVariables(void)
{
	// FIXME
        locVmax = 0.0; 
	totVmax = 0.0;
        locHalosSize[0] = 0; locHalosSize[1] = 0;
	nLocParts[0] = 0; nLocParts[1] = 0;

	// We will do a pairwise comparison of the catalogs
	locParts.resize(2);
	locHalos.resize(2);

	GlobalGrid[0].Init(nGrid, boxSize);
	GlobalGrid[1].Init(nGrid, boxSize);

	sizeHalo = sizeof(Halo);
	sizePart = sizeof(unsigned long long int);
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


void CleanMemory(int iCat)
{
	locHalos[iCat].clear();
	locHalos[iCat].shrink_to_fit();

		for (int iH = 0; iH < nLocHalos[iCat]; iH++)
		{
			for (int iT = 0; iT < 6; iT++)
			{
				locParts[iCat][iH][iT].clear();
				locParts[iCat][iH][iT].shrink_to_fit();
			}
				
			locParts[iCat][iH].clear();
			locParts[iCat][iH].shrink_to_fit();
		}

		locParts[iCat].clear();
		locParts[iCat].shrink_to_fit();
};


void ShiftHalosAndParts()
{
	if (locTask == 0)
		cout << "Shifting halos and particles..." << endl;

	CleanMemory(0);

	nLocHalos[0] = nLocHalos[1];
	nLocParts[0] = nLocParts[1];

	locHalos[0] = locHalos[1];
	locParts[0].resize(nLocHalos[0]);

	for (int iH = 0; iH < nLocHalos[0]; iH++)
	{
		locParts[0][iH].resize(nPTypes);

		for (int iT = 0; iT < 6; iT++)
			locParts[0][iH][iT] = locParts[1][iH][iT];
	}
};


float VectorModule(float *V)
{
	return sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
};

