#include <mpi.h>

#include <iostream>
#include <string>

#include "IOSettings.h"
#include "Halo.h"
#include "Particle.h"
#include "general.h"

using namespace std;



int main(int argv, char **argc)
{
	int nCpus = 0;	// Number of cpus used for the analysis - each file (part or halo) is split into nCpus parts
	int nSnap = 0;	// Total number of snapshots - e.g. from 0 to 54 
	int iStep = 0;	// Step of the iteration	
	
	int iHalo = 0, jHalo = 0, kHalo = 0;
	
	size_t locHaloSize = 0;
	size_t sizeP = 0;

	string partFilesInfo;
	string haloFilesInfo;

	IOSettings SettingsIO;

	MPI_Init(&argv, &argc);
	MPI_Comm_rank(MPI_COMM_WORLD, &locTask);
 	MPI_Comm_size(MPI_COMM_WORLD, &totTask);

	locHalosSize = 0;
	totHalosSize = 0;

	nLocHalos = 10;

	// Each task could read more than one file, this ensures it only reads adjacent snapshots
	SettingsIO.DistributeFilesAmongTasks();

	SettingsIO.ReadHalos();

	locHalos = new Halo[nLocHalos];
	
	Particle *p;
	p = new Particle[10];

	locHalos[0].ID = 1 + locTask;
	locHalos[1].Mtot = 10 * (locTask + 1);

	for (iHalo = 0; iHalo < nLocHalos; iHalo++)
	{
		locHalos[iHalo].Part = new Particle[iHalo+1];

		locHalosSize = sizeof(locHalos[iHalo]);
		totHalosSize += locHalosSize;

		sizeP = sizeof(Particle[0]);

		cout << "ThisTask: " << locTask << ", halo=" << iHalo << ", HaloSize: " << locHalosSize
			<< ", totSize : " << totHalosSize << ", pSize: " << sizeP << endl; 
	}

	// Read the halo files - one per task

	// Read the particle files - one per task

	// Get the maximum halo velocity to compute buffer size

	// Split the files among the tasks:
	// - compute the optimal size of the splitting region according to the buffer size
	// - if ntask > nfiles read in then distribute the halos 
	
	MPI_Finalize();
}


