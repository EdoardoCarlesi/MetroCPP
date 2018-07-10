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

	// TODO Make these parameters readable from an input file
	int totCat = 2;	
	boxSize = 1.0e+5; 	// Using kpc/h units
	nGrid = 100;	
	nChunksPerFile = 1;
	nPTypes = 6;
	GeneralMethods.dMaxFactor = 1.5;

	MPI_Init(&argv, &argc);
	MPI_Comm_rank(MPI_COMM_WORLD, &locTask);
 	MPI_Comm_size(MPI_COMM_WORLD, &totTask);
	
	InitLocVariables();

	// TODO make this readable from input file
	SettingsIO.pathInput = "/home/eduardo/CLUES/DATA/FullBox/01/";
	SettingsIO.catFormat = "AHF_halos";
	SettingsIO.thisPath = "/home/eduardo/CLUES/MetroC++/";
	SettingsIO.haloSuffix = "AHF_halos";
	SettingsIO.partSuffix = "AHF_particles";
	SettingsIO.haloPrefix = "snapshot_";
	SettingsIO.partPrefix = "snapshot_";

	SettingsIO.Init();

	// We are assuming that each task reads more than one file
	SettingsIO.DistributeFilesAmongTasks();

	iUseCat = 0;
	SettingsIO.ReadHalos();
	SettingsIO.ReadParticles();	

	// TODO some load balancing
	// Split the files among the tasks:
	// - if ntask > nfiles read in then distribute the halos 

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

		//MPI_Barrier(MPI_COMM_WORLD);

		// After reading in the second halo catalog, each task finds out which buffer nodes it needs to request to the other tasks
		// The nodes are located on grid 1 based on the distribution of the nodes on grid 0
		GlobalGrid[1].FindBufferNodes(GlobalGrid[0].locNodes);	// TODO Check with globalGrid[1] there might be some allocation
									// problem, or initialization!!!!!
		// Communicate the list of the buffer nodes to be exchanged among every task
		//CommTasks.ExchangeBuffers();

		// Now exchange the halos in the requested buffer zones among the different tasks
		//CommTasks.BufferSendRecv();

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
	
#ifdef TEST
#endif
		// Now shift the halo catalog from 1 to 0, and clean the buffers
		ShiftHalosPartsGrids();
		CleanMemory(1);
	}
	
	if (locTask == 0)
		cout << "The loop on halo and particle catalogs has finished." << endl;

	// Retrieve some informations on the grid - sanity check
	//GlobalGrid.Info();

	// Get the maximum halo velocity to compute buffer size


	//cout << "Finished on task=" << locTask << endl;

	CleanMemory(0);

#ifdef TEST_BLOCK
#endif
	MPI_Finalize();
	
	if (locTask == 0)	
		cout << "MPI finalized, memory cleaned. Exiting the program." << endl;

	//exit(0);
}


