#include <mpi.h>

#include "Comm.h"

// MPI communication routines that take care of broadcasting, swapping, packing and send/recv of messages 

/* This function computes the size of all the packages that need to be communicated across the tasks 
 * Halos within a "buffer" volume are MPI_Pack-ed and sent to the requesting task
 */
void Comm::ComputeBufferSize()
{
};


void Comm::BufferSendRecv()
{
};


//void MpiComm::

