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

Grid GlobalGrid;

vector<vector<Halo>> locHalos;
vector <Halo> tmpHalos;

size_t locHalosSize[2];
size_t tmpHalosSize;

int iNumCat; // Halo catalog number in use, from 0 to N
int iUseCat; // Refers to 0 or 1 depending on the snapshot being used
int nTypePart;
int nTotHalos[2];
int nLocHalos[2];

vector<vector<vector<vector<unsigned long long int>>>> locParts;
vector<vector<unsigned long long int>> tmpParts;	

size_t locPartsSize[2];
size_t tmpPartsSize;

size_t sizeHalo;
size_t sizePart;

int nPTypes;
int nTmpParts;
int nLocParts[2];

float totVmax;
float locVmax;
float bufferThickness;
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
        tmpHalosSize = 0; 

	// We will do a pairwise comparison of the catalogs
	locParts.resize(2);
	locHalos.resize(2);

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


void CleanMemory()
{
	locHalos[iUseCat].clear();
	locHalos[iUseCat].shrink_to_fit();

		for (int iH = 0; iH < nLocHalos[iUseCat]; iH++)
		{
			for (int iT = 0; iT < 6; iT++)
			{
				locParts[iUseCat][iH][iT].clear();
				locParts[iUseCat][iH][iT].shrink_to_fit();
			}
				
			locParts[iUseCat][iH].clear();
			locParts[iUseCat][iH].shrink_to_fit();
		}

		locParts[iUseCat].clear();
		locParts[iUseCat].shrink_to_fit();

};


float VectorModule(float *V)
{
	return sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
};

