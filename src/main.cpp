#include <mpi.h>

#include <vector>
#include <iostream>
#include <string>
#include <ctime>

#include "general.h"
#include "Communication.h"
#include "IOSettings.h"
#include "Halo.h"
#include "Methods.h"

using namespace std;



int main(int argv, char **argc)
{
	IOSettings SettingsIO;
	Communication CommTasks;
	Methods GeneralMethods;

	// Make these parameters readable from an input file
	int totCat = 3;	
	float boxSize = 1.0e+5; 	// Using kpc/h units
	int nGrid = 100;	
	nChunksPerFile = 2;
	nPTypes = 6;

	MPI_Init(&argv, &argc);
	MPI_Comm_rank(MPI_COMM_WORLD, &locTask);
 	MPI_Comm_size(MPI_COMM_WORLD, &totTask);
	
	InitLocVariables();

	GlobalGrid[0].Init(nGrid, boxSize);
	GlobalGrid[1].Init(nGrid, boxSize);

	// TODO make this readable from input file
	SettingsIO.pathInput = "/home/eduardo/CLUES/DATA/FullBox/01/";
	SettingsIO.catFormat = "AHF_halos";
	SettingsIO.thisPath = "/home/eduardo/CLUES/MetroC++/";
	SettingsIO.haloSuffix = "AHF_halos";
	SettingsIO.partSuffix = "AHF_particles";
	SettingsIO.haloPrefix = "snapshot_";
	SettingsIO.partPrefix = "snapshot_";

	SettingsIO.Init();

	// Each task could read more than one file, this ensures it only reads adjacent snapshots
	SettingsIO.DistributeFilesAmongTasks();
	GeneralMethods.dMaxFactor = 1.5;

	iUseCat = 0;
	SettingsIO.ReadHalos();
	SettingsIO.ReadParticles();	

	// Now every task knows which nodes belongs to which task
	CommTasks.BroadcastAndGatherGrid();

	for (iNumCat = 1; iNumCat < totCat; iNumCat++)
	{
		clock_t iniTime = clock();

		iUseCat = 1;
		SettingsIO.ReadHalos();
		SettingsIO.ReadParticles();	

		// Now exchange the halos in the requested buffer zones among the different tasks
		//CommTasks.BufferSendRecv();

		// Now every task knows which nodes belongs to which task
		CommTasks.BroadcastAndGatherGrid();

		if (locTask == 0)
			cout << "Finding halo progentors, forwards..." << flush ;
	
		GeneralMethods.FindProgenitors(0, 1);

		if (locTask == 0)
			cout << "\nFinding halo progentors, backwards..." << flush ;
	
		GeneralMethods.FindProgenitors(1, 0);

		clock_t endTime = clock();
		double elapsed = double(endTime - iniTime) / CLOCKS_PER_SEC;

		if (locTask == 0)
			cout << "\nDone in " << elapsed << "s. " << endl;
	
		// Now shift the halo catalog from 1 to 0, and clean the buffers
		ShiftHalosAndParts();
		CleanMemory(1);

#ifdef TEST_BLOCK
#endif
	}
	
	if (locTask == 0)
		cout << "The loop on halo and particle catalogs has finished." << endl;

	// Retrieve some informations on the grid - sanity check
	//GlobalGrid.Info();

	// Get the maximum halo velocity to compute buffer size

	// Split the files among the tasks:
	// - compute the optimal size of the splitting region according to the buffer size
	// - if ntask > nfiles read in then distribute the halos 

	//cout << "Finished on task=" << locTask << endl;

	CleanMemory(0);

	int end = MPI_Finalize();
	
	if (locTask == 0)	
		cout << "MPI finalized, memory cleaned. Exiting the program." << endl;

	exit(0);
}


