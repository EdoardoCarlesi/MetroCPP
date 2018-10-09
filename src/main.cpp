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

	InitTrees(nSnaps);

	string strRunMode;

	if (runMode == 0)
		strRunMode = " ---> Merger tree computation only.\n";
	else if (runMode == 1)	
		strRunMode = " ---> Post processing the trees only.\n";
	else if (runMode == 2)
		strRunMode = " ---> Merger tree computation and post processing.\n";
	else {
		if (locTask == 0)
		{
			cout << "RunMode unknown. Choose 0, 1 or 2." << endl;
			exit(0);
		}
	};

	if (locTask == 0)
		cout << "Running the code in mode: " << runMode << strRunMode << endl;

	/* If running in MTree only or MTree + Postprocessing */
	if (runMode == 0 || runMode == 2)
	{
		/* Overrides the config file settings */
		nTreeChunks = totTask;
	
		/* We are assuming that each task reads more than one file. TODO load balancing */
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

		if (locTask == 0)
			cout << "Starting to loop on " << nSnaps << " halo and particle files." << endl;

		/* Loop on halo and particle catalogs */
		for (iNumCat = 1; iNumCat < nSnaps; iNumCat++)
		{
			clock_t iniTime = clock();
		
			iUseCat = 1;
			SettingsIO.ReadHalos();
			SettingsIO.ReadParticles();	
#ifndef ZOOM
			/* Now every task knows which nodes belongs to which task */
			CommTasks.BroadcastAndGatherGrid();

			/* After reading in the second halo catalog, each task finds out which nodes it gets from the other tasks
			 * The nodes are located on grid 1 based on the distribution of the nodes on grid 0 */
			GlobalGrid[1].FindBufferNodes(GlobalGrid[0].locNodes);	
#endif
			/* Now exchange the halos in the requested buffer zones among the different tasks */
			CommTasks.BufferSendRecv();

			if (locTask == 0)
				cout << "Finding halo progentors, forwards..." << flush ;
		
			/* This function also allocates the MergerTrees */
			FindProgenitors(0, 1);

			clock_t endTime = clock();
			double elapsed = double(endTime - iniTime) / CLOCKS_PER_SEC;

			if (locTask == 0)
				cout << "\nDone in " << elapsed << "s. " << endl;
		
			MPI_Barrier(MPI_COMM_WORLD);
			CommTasks.SyncOrphanHalos();

			if (locTask == 0)
				cout << "\nFinding halo progentors, backwards..." << flush ;

			iniTime = clock();
	
			FindProgenitors(1, 0);

			endTime = clock();
			elapsed = double(endTime - iniTime) / CLOCKS_PER_SEC;
	
			if (locTask == 0)
				cout << "\nDone in " << elapsed << "s. " << endl;

			CleanTrees(iNumCat);

			/* Now shift the halo catalog from 1 to 0, and clean the buffers */
			ShiftHalosPartsGrids();

			CleanMemory(1);
		}	/* The trees have now been built */

		//DebugTrees();

		if (locTask == 0)
			cout << "The loop on halo and particle catalogs has finished." << endl;

		MPI_Barrier(MPI_COMM_WORLD);
			SettingsIO.WriteTrees();
	
		CleanMemory(0);

	}	/* If running the tree and / or post processing mode only */
	
	/* Load in trees & halo catalogs */
	if (runMode == 1)
	{
		iNumCat = 0;
		iUseCat = 0;
		SettingsIO.ReadHalos();

		for (iNumCat = 1; iNumCat < nSnaps; iNumCat++)
		{
			SettingsIO.ReadHalos();
			SettingsIO.ReadTrees();

			// Do something to add orphan trees to the next step
		}
	}

	/* Proceed with smoothing & interpolating the MAH of single halos */
	if (runMode == 1 || runMode == 2)
	{
		/* 
			POST PROCESSING STUFF:
				- interpolate lost halo masses 
				- compute local gravitational field to find missing halo's positions
				- smooth over the mass: include (exclude) transient halos, flybys, subhalos outside Rvir
		*/
	}

	MPI_Finalize();
	
	if (locTask == 0)	
		cout << "MPI finalized, memory cleaned. Exiting the program." << endl;

	exit(0);
}


