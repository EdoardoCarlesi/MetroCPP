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
#include "Cosmology.h"

using namespace std;



int main(int argv, char **argc)
{
	IOSettings SettingsIO;
	Communication CommTasks;
	Cosmology Cosmo;

	string configFile = argc[1];

	/* Read configuration file and initialize variables */
	SettingsIO.ReadConfigFile(configFile);

	/* Fire off MPI */
	MPI_Init(&argv, &argc);
	MPI_Comm_rank(MPI_COMM_WORLD, &locTask);
 	MPI_Comm_size(MPI_COMM_WORLD, &totTask);
                
	/* Error handler */
	//MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN); 
	//MPI::COMM_WORLD.Set_errhandler ( MPI::ERRORS_THROW_EXCEPTIONS );
	//MPI::COMM_WORLD.Set_errhandler ( MPI::ERRORS_RETURN);

	/* These are the local variables and I/O settings that are defined for each task*/
	InitLocVariables();

	SettingsIO.Init();

	InitTrees(nSnapsUse);
	SettingsIO.SetCosmology(&Cosmo);

	string strRunMode;

	/* We are assuming that each task reads more than one file. TODO load balancing */
	SettingsIO.DistributeFilesAmongTasks();

	if (locTask == 0)
	{
		cout << endl;
		cout << "\t\t=========================================" << endl;
#ifdef ZOOM
		cout << "\t\t=========== ZOOM OPERATION MODE =========" << endl;
#else
		cout << "\t\t========= FULL BOX OPERATION MODE =======" << endl;
#endif
		cout << "\t\t=========================================" << endl;
		cout << endl;
	}

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
			cout << "Starting to loop on " << nSnapsUse << " halo and particle files." << endl;

		/* Loop on halo and particle catalogs */
		for (iNumCat = 1; iNumCat < nSnapsUse; iNumCat++)
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

			/* Now exchange the halos in the requested buffer zones among the different tasks.
			 * In zoom mode we send ALL halos to ALL tasks */
			//CommTasks.BufferSendRecv();

			if (locTask == 0)
				cout << "Finding halo progentors, forwards..." << flush ;
		
			/* Forward halo connections
			 * This function also allocates the MergerTrees */
			FindProgenitors(0, 1);
			MPI_Barrier(MPI_COMM_WORLD);

			clock_t endTime = clock();
			double elapsed = double(endTime - iniTime) / CLOCKS_PER_SEC;

			if (locTask == 0)
				cout << "\nDone in " << elapsed << "s. " << endl;
		
#ifdef ZOOM
			/* Orphan halo candidates need to be communicated in zoom mode only. 
			 * In fullbox mode they are taken care of in the FindProgenitors function */
			CommTasks.SyncOrphanHalos();
#endif

			if (locTask == 0)
				cout << "\nFinding halo progentors, backwards..." << flush ;

			iniTime = clock();
	
			/* Backward halo connections */
			FindProgenitors(1, 0);
			MPI_Barrier(MPI_COMM_WORLD);
	
			endTime = clock();
			elapsed = double(endTime - iniTime) / CLOCKS_PER_SEC;
	
			if (locTask == 0)
				cout << "\nDone in " << elapsed << "s. " << endl;

			CleanTrees(iNumCat);

			/* Now shift the halo catalog from 1 to 0, and clean the buffers */
#ifndef ZOOM
			CommTasks.CleanBuffer();
#endif
			ShiftHalosPartsGrids();
			SettingsIO.WriteTree(iNumCat); 	// TODO write trees at each step
			
		}	/* Finish: the trees have now been built for this step */

		//SettingsIO.WriteTrees(); 	// TODO write trees at each step
		if (locTask == 0)
			cout << "The loop on halo and particle catalogs has finished." << endl;
	
		CleanMemory(0);

	}	/* If running the tree and / or post processing mode only */
	
	/* Load in trees & halo catalogs */
	if (runMode == 1)
	{
		iNumCat = 0;	iUseCat = 0;

		SettingsIO.ReadHalos();

#ifdef ZOOM	
		/* We need to scatter catalogs through the tasks only when in ZOOM mode	*/
		CommTasks.BufferSendRecv();
#else		
		/* Communicate grid info across all tasks */
		CommTasks.BroadcastAndGatherGrid();
#endif
		for (iNumCat = 1; iNumCat < nSnapsUse; iNumCat++)
		{
			/* Initialize descendant halos */
			iUseCat = 0;			
			CommTasks.SyncIndex();
			SettingsIO.ReadTrees();
			AssignDescendant();

			// TODO: At this point we could fix halo masses and position using some interpolation scheme
			
			iUseCat = 1;
			SettingsIO.ReadHalos();
#ifdef ZOOM
			CommTasks.BufferSendRecv();
			CommTasks.SyncOrphanHalos();
			AssignProgenitor();
#else
			CommTasks.BroadcastAndGatherGrid();
			GlobalGrid[1].FindBufferNodes(GlobalGrid[0].locNodes);	
			CommTasks.BufferSendRecv();
			CommTasks.SyncIndex();
			AssignProgenitor();
			CommTasks.CleanBuffer();
#endif

			ShiftHalosPartsGrids();
		}

		SettingsIO.WriteTrees();

		// TODO: once the MAHs have been computed, we need to smooth over 
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


