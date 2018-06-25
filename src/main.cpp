#include <mpi.h>

#include <iostream>
#include <string>

#include "Communication.h"
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
	
	float boxSize = 1.0e+5; int nGrid = 100;	// Unsing kpc/h units

	string partFilesInfo;
	string haloFilesInfo;

	IOSettings SettingsIO;
	Communication CommTasks;

	nChunksPerFile = 8;

	MPI_Init(&argv, &argc);
	MPI_Comm_rank(MPI_COMM_WORLD, &locTask);
 	MPI_Comm_size(MPI_COMM_WORLD, &totTask);
	
	InitLocVariables();
	GlobalGrid.Init(nGrid, boxSize);

	// Each task could read more than one file, this ensures it only reads adjacent snapshots
	SettingsIO.DistributeFilesAmongTasks();
	SettingsIO.pathInput = "/home/eduardo/CLUES/DATA/FullBox/01/";

	string fileRoot = "snapshot_054.000";
	string fileSuffHalo = ".z0.000.AHF_halos";
	string fileSuffPart = ".z0.000.AHF_particles";
	string nameTask = to_string(locTask);

	SettingsIO.urlTestFileHalo = fileRoot + nameTask + fileSuffHalo;
	SettingsIO.urlTestFilePart = fileRoot + nameTask + fileSuffPart;

	// Read the halo files - one per task
	// This is also assigning the halo list to each grid node
	SettingsIO.ReadHalos();

	// Determine which task holds which sub-volume of the box
	//GlobalGrid.FindPatchOnTask();

	// Now every task knows which nodes belongs to which task
	CommTasks.BroadcastAndGatherGrid();

	//GlobalGrid.Info();

	// Read the particle files - one per task
	//SettingsIO.ReadParticles();	

	// Now exchange the halos in the requested buffer zones among the different tasks
	//CommTasks.BufferSendRecv();

	// Get the maximum halo velocity to compute buffer size

	// Split the files among the tasks:
	// - compute the optimal size of the splitting region according to the buffer size
	// - if ntask > nfiles read in then distribute the halos 

	//cout << "Finished on task=" << locTask << endl;

#ifdef TEST_BLOCK
	locHalos.clear();
	locHalos.shrink_to_fit();
#endif

	MPI_Finalize();
	
	if (locTask == 0)	
		cout << "Done." << endl;

	exit(0);
}


