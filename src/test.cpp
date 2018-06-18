#include <mpi.h>
#include <iostream>

#include "IOSettings.h"
#include "Halo.h"
#include "Particle.h"
#include "general.h"

using namespace std;

int main(int argv, char **argc)
{
	MPI_Init(&argv, &argc);
	MPI_Comm_rank(MPI_COMM_WORLD, &locTask);
 	MPI_Comm_size(MPI_COMM_WORLD, &totTask);

	Halo *locHalo;

	locHalo = new Halo[2];

	locHalo[0].Info();

	// Read the halo files - one per task

	// Read the particle files - one per task

	// Get the maximum halo velocity to compute buffer size

	// Split the files among the tasks:
	// - compute the optimal size of the splitting region according to the buffer size
	// - if ntask > nfiles read in then distribute the halos 


	// Remove the buffer halos AND the halos at step n-1, load the halos at step n+1 AND re-compute the buffer
	
	MPI_Finalize();
}


