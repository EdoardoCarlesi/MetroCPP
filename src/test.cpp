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
	int nCpus = 0;	// Number of cpus used for the analysis - each file (part or halo) is split into nCpus parts
	int nSnap = 0;	// Total number of snapshots - e.g. from 0 to 54 
	int iStep = 0;	// Step of the iteration	
	
	int iHalo = 0, jHalo = 0, kHalo = 0;
		
	vector <int> nCommon;
	size_t locHaloSize = 0;
	size_t sizeP = 0;
	
	float boxSize = 1.0e+5; int nGrid = 100;	// Unsing kpc/h units

	string partFilesInfo;
	string haloFilesInfo;

	IOSettings SettingsIO;
	Communication CommTasks;
	Methods GeneralMethods;

	nChunksPerFile = 8;

	MPI_Init(&argv, &argc);
	MPI_Comm_rank(MPI_COMM_WORLD, &locTask);
 	MPI_Comm_size(MPI_COMM_WORLD, &totTask);
	
	InitLocVariables();
	GlobalGrid.Init(nGrid, boxSize);

	// Each task could read more than one file, this ensures it only reads adjacent snapshots
	//SettingsIO.DistributeFilesAmongTasks();

	// TODO make this readable from input file
	SettingsIO.pathInput = "/home/eduardo/CLUES/DATA/FullBox/01/";
	SettingsIO.catFormat = "AHF_halos";
	SettingsIO.thisPath = "/home/eduardo/CLUES/MetroC++/";
	SettingsIO.haloSuffix = "AHF_halos";
	SettingsIO.partSuffix = "AHF_particles";
	SettingsIO.haloPrefix = "snapshot_";
	SettingsIO.partPrefix = "snapshot_";

	string fileRoot = "snapshot_054.000";
	string fileSuffHalo = ".z0.000.AHF_halos";
	string fileSuffPart = ".z0.000.AHF_particles";
	string nameTask = to_string(locTask);

	SettingsIO.urlTestFileHalo = fileRoot + nameTask + fileSuffHalo;
	SettingsIO.urlTestFilePart = fileRoot + nameTask + fileSuffPart;

	SettingsIO.Init();

	SettingsIO.DistributeFilesAmongTasks();

	int totCat = 1;	// Only read the 

	for (int iCat = 0; iCat < totCat; iCat++)
	{
		// Read the halo files - one per task
		// This is also assigning the halo list to each grid node
		SettingsIO.ReadHalos();

		// Read the particle files - one per task
		SettingsIO.ReadParticles();	

		nCommon.resize(locHalos[0].nTypes);
		//cout << "Print part= " << locParts[0][1][1] << endl;

		clock_t iniTime = clock();
		nCommon = GeneralMethods.CommonParticles(locParts[0], locParts[0]);
		clock_t endTime = clock();

		double elapsed = double(endTime - iniTime) / CLOCKS_PER_SEC;
	
		//cout << locParts[0][0].first << " " << locParts[0][0].second << " on task " << locTask << endl; 
		//cout << "Npart: " << locParts[0].size() << " on task " << locTask << endl; 
		//cout << "Npart: " << locHalos[0].nPart << " on task " << locTask << endl; 

		cout << "Process took " << elapsed << " s on task "<< locTask << endl;
	}


	// Now every task knows which nodes belongs to which task
	//CommTasks.BroadcastAndGatherGrid();

	// Retrieve some informations on the grid - sanity check
	//GlobalGrid.Info();


	// Now exchange the halos in the requested buffer zones among the different tasks
	//CommTasks.BufferSendRecv();

	// Get the maximum halo velocity to compute buffer size

	// Split the files among the tasks:
	// - compute the optimal size of the splitting region according to the buffer size
	// - if ntask > nfiles read in then distribute the halos 

	//cout << "Finished on task=" << locTask << endl;

	CleanMemory();

	MPI_Finalize();
	
	if (locTask == 0)	
		cout << "Done." << endl;

	exit(0);
}


