#include <mpi.h>

#include <vector>
#include <iostream>
#include <string>
#include <ctime>

#include "global_vars.h"
#include "utils.h"

#include "Communication.h"
#include "IOSettings.h"
#include "Halo.h"

using namespace std;



int main(int argv, char **argc)
{
	IOSettings SettingsIO;
	Communication CommTasks;

	// TODO Make these parameters readable from an input file
	int totCat = 2;	
	boxSize = 1.0e+5; 	// Using kpc/h units
	nGrid = 100;	
	nChunksPerFile = 2;
	dMaxFactor = 1.5;

	string configFile = argc[1];

	MPI_Init(&argv, &argc);
	MPI_Comm_rank(MPI_COMM_WORLD, &locTask);
 	MPI_Comm_size(MPI_COMM_WORLD, &totTask);

	InitLocVariables();
	//SettingsIO.ReadConfigFile(configFile);

	// TODO make this readable from input file
	SettingsIO.pathMetroCpp = "/home/eduardo/CLUES/MetroC++/";
	SettingsIO.pathInput = "/home/eduardo/CLUES/DATA/FullBox/01/";
	SettingsIO.catFormat = "AHF_halos";
	SettingsIO.haloSuffix = "AHF_halos";
	SettingsIO.partSuffix = "AHF_particles";
	SettingsIO.haloPrefix = "snapshot_";
	SettingsIO.partPrefix = "snapshot_";

	SettingsIO.Init();

	// We are assuming that each task reads more than one file
	SettingsIO.DistributeFilesAmongTasks();

	/* This is a global variable */
	iUseCat = 0;

	/* Read particles and catalogs */
	SettingsIO.ReadHalos();
	SettingsIO.ReadParticles();	

	// TODO some load balancing
	// Split the files among the tasks in case ntask > nfiles to be read in

	/* Now every task knows which subvolumes of the box belong to which task */
	CommTasks.BroadcastAndGatherGrid();

	/* Loop on halo and particle catalogs */
	for (iNumCat = 1; iNumCat < totCat; iNumCat++)
	{
		clock_t iniTime = clock();

		iUseCat = 1;
		SettingsIO.ReadHalos();
		SettingsIO.ReadParticles();	

		// Now every task knows which nodes belongs to which task
		CommTasks.BroadcastAndGatherGrid();

		// After reading in the second halo catalog, each task finds out which buffer nodes it needs to request to the other tasks
		// The nodes are located on grid 1 based on the distribution of the nodes on grid 0
		GlobalGrid[1].FindBufferNodes(GlobalGrid[0].locNodes);	

		// Now exchange the halos in the requested buffer zones among the different tasks
		CommTasks.BufferSendRecv();

		MPI_Barrier(MPI_COMM_WORLD);

		if (locTask == 0)
			cout << "Finding halo progentors, forwards..." << flush ;
	
		FindProgenitors(0, 1);

		clock_t endTime = clock();
		double elapsed = double(endTime - iniTime) / CLOCKS_PER_SEC;

		if (locTask == 0)
			cout << "\nDone in " << elapsed << "s. " << endl;
	

#ifdef TEST
		clock_t iniTime = clock();

		if (locTask == 0)
			cout << "\nFinding halo progentors, backwards..." << flush ;

	
		FindProgenitors(1, 0);

		clock_t endTime = clock();
		double elapsed = double(endTime - iniTime) / CLOCKS_PER_SEC;

		if (locTask == 0)
			cout << "\nDone in " << elapsed << "s. " << endl;
	
#endif
		// Now shift the halo catalog from 1 to 0, and clean the buffers
		ShiftHalosPartsGrids();
		CleanMemory(1);
	}
	
#ifdef TEST_BLOCK
	if (locTask == 0)
		cout << "The loop on halo and particle catalogs has finished." << endl;

	// Retrieve some informations on the grid - sanity check
	//GlobalGrid.Info();

	// Get the maximum halo velocity to compute buffer size


	//cout << "Finished on task=" << locTask << endl;

#endif
	CleanMemory(0);

	MPI_Finalize();
	
	if (locTask == 0)	
		cout << "MPI finalized, memory cleaned. Exiting the program." << endl;

	//exit(0);
}


