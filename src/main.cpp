/* The main file is a wrapper for all the functions that need to be invoked to 
 * compute the merger trees. 
 * First, we read an input file containing the necessary settings.
 * Then halo/particle files are distributed and read in by the tasks.
 * Buffers along the borders of the subvolume contained by each task are then computed and distributed across the tasks.
 * Halo merger trees are then computed consistently on each individual task. */

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
#include "MergerTree.h"

using namespace std;



int main(int argv, char **argc)
{
	IOSettings SettingsIO;
	Communication CommTasks;

	string configFile = argc[1];

	/* Read configuration file and initialize variables */
	SettingsIO.ReadConfigFile(configFile);

	/* Fire off MPI */
	MPI_Init(&argv, &argc);
	MPI_Comm_rank(MPI_COMM_WORLD, &locTask);
 	MPI_Comm_size(MPI_COMM_WORLD, &totTask);

	/* These are the local variables and I/O settings that are defined for each task*/
	InitLocVariables();
	SettingsIO.Init();

	// We are assuming that each task reads more than one file. TODO load balancing
	SettingsIO.DistributeFilesAmongTasks();

	/* This is a global variable */
	iUseCat = 0;

	/* Read particles and catalogs */
#ifdef ZOOM
	if (locTask ==0)
	{
#endif
		SettingsIO.ReadHalos();
		SettingsIO.ReadParticles();	
#ifdef ZOOM
	}
#endif

#ifdef VERBOSE
	if (locTask == 0)
		SettingsIO.CheckStatus();
#endif

#ifdef ZOOM
	/* Broadcast all the halos on all tasks */
	CommTasks.BufferSendRecv();
#else
	/* Now every task knows which subvolumes of the box belong to which task */
	CommTasks.BroadcastAndGatherGrid();
#endif

	int nUseCat = 10;	// THIS IS A LOCAL VARIABLE used for TEST only

	// TODO put all these functions somewhere into the MergerTree file
	locCleanTrees.resize(nUseCat-1);
	allHalos.resize(nUseCat);
	copy(locHalos[0].begin(), locHalos[0].end(), back_inserter(allHalos[iNumCat]));

	/* Loop on halo and particle catalogs */
	for (iNumCat = 1; iNumCat < nUseCat; iNumCat++)
	{
		clock_t iniTime = clock();

		iUseCat = 1;
		SettingsIO.ReadHalos();
		SettingsIO.ReadParticles();	

#ifndef ZOOM
		/* Now every task knows which nodes belongs to which task */
		CommTasks.BroadcastAndGatherGrid();

		/* After reading in the second halo catalog, each task finds out which nodes it gets from the other tasks
		   The nodes are located on grid 1 based on the distribution of the nodes on grid 0 */
		GlobalGrid[1].FindBufferNodes(GlobalGrid[0].locNodes);	
#endif

		/* Now exchange the halos in the requested buffer zones among the different tasks */
		CommTasks.BufferSendRecv();

		//MPI_Barrier(MPI_COMM_WORLD);

		if (locTask == 0)
			cout << "Finding halo progentors, forwards..." << flush ;
		
		/* This function also allocates the MergerTrees */
		FindProgenitors(0, 1);

		if (locTask == 0)
			cout << "\nFinding halo progentors, backwards..." << flush ;
	
		FindProgenitors(1, 0);

		clock_t endTime = clock();
		double elapsed = double(endTime - iniTime) / CLOCKS_PER_SEC;

		if (locTask == 0)
			cout << "\nDone in " << elapsed << "s. " << endl;
	
		CleanTrees(iNumCat);

		// Now shift the halo catalog from 1 to 0, and clean the buffers
		ShiftHalosPartsGrids();

		// TODO move this copy somewhere else
		copy(locHalos[1].begin(), locHalos[1].end(), back_inserter(allHalos[iNumCat]));
		CleanMemory(1);
	}
	
	if (locTask == 0)
		cout << "The loop on halo and particle catalogs has finished." << endl;


	// The trees have now been built


#ifndef ZOOM
	// Retrieve some informations on the grid - sanity check
	//GlobalGrid.Info();
#endif

	//MPI_Barrier(MPI_COMM_WORLD);

	DebugTrees();

	CleanMemory(0);

	MPI_Finalize();
	
	if (locTask == 0)	
		cout << "MPI finalized, memory cleaned. Exiting the program." << endl;

	exit(0);
}


