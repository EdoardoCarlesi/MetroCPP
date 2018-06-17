#include <mpi.h>

#include "Halo.h"
#include "Particle.h"
#include "general.h"

int main()
{

	MPI::Init();
	LocTask = MPI::COMM_WORLD.Get_rank();
	TotTask = MPI::COMM_WORLD.Get_size();
	MPI::Finalize();

}


