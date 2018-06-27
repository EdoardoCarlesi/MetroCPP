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

vector <Halo> locHalos;
size_t locHalosSize;
size_t totHalosSize;

int nTypePart;
int nTotHalos;
int nLocHalos;
int iLocHalos;

//vector <vector <Particle>> locParts;
//vector <unordered_map<unsigned long long int, short int>> locParts;
vector <vector<vector<unsigned long long int>>> locParts;	// numbers of particles are stored by particle type
size_t locPartsSize;
size_t totPartsSize;

size_t sizeHalo;
size_t sizePart;

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
	sizePart = sizeof(unsigned long long int) + sizeof(short int);	

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


void CleanMemory()
{

	locHalos.clear();
	locHalos.shrink_to_fit();

		for (int iH = 0; iH < nLocHalos; iH++)
		{
			for (int iT = 0; iT < 6; iT++)
			{
				locParts[iH][iT].clear();
				locParts[iH][iT].shrink_to_fit();
			}
				
			locParts[iH].clear();
			locParts[iH].shrink_to_fit();
		}

		locParts.clear();
		locParts.shrink_to_fit();

};


float VectorModule(float *V)
{
	return sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
};

