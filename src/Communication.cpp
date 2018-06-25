#include <mpi.h>
#include <iostream>
#include <vector>

#include "general.h"
#include "Grid.h"
#include "Halo.h"
#include "Communication.h"

using namespace std;


// MPI communication routines that take care of broadcasting, swapping, packing and send/recv of messages 

/* This function computes the size of all the packages that need to be communicated across the tasks 
 * Halos within a "buffer" volume are MPI_Pack-ed and sent to the requesting task
 */
void Communication::ComputeBufferSize()
{
};


/*
void Comm::Optimize()	// Write some function that optimizes the memory distribution among the tasks
{};
*/

void Communication::BroadcastAndGatherGrid()
{
	size_t gridSize = GlobalGrid.nNodes;
	vector <vector <int>> tmpTaskOnGridNode;
	//vector <int> tmpTaskOnGridNode;
 
	//if (locTask == 0)
	{
		tmpTaskOnGridNode.resize(totTask);

		//tmpTaskOnGridNode.resize(gridSize);
		for (int i = 0; i < totTask; i++)
			tmpTaskOnGridNode[i].resize(gridSize);
	}

//MPI_Reduce(&GlobalGrid.taskOnGridNode[0], &tmpTaskOnGridNode[locTask][0], gridSize, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
//MPI_Allreduce(&GlobalGrid.taskOnGridNode[0], &tmpTaskOnGridNode[0], gridSize, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//MPI_Allreduce(&GlobalGrid.taskOnGridNode[0], &tmpTaskOnGridNode[0], gridSize, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//MPI_Gather(&GlobalGrid.taskOnGridNode[0], gridSize, MPI_INT, &tmpTaskOnGridNode[locTask][0], gridSize, MPI_INT, 0, MPI_COMM_WORLD);

#ifdef TEST_BLOCK
	if (locTask == 0)
	{
		int *iTasksAll; iTasksAll = new int[totTask];
		//GlobalGrid.globalTaskOnGridNode.resize(gridSize);
		
		for (int j = 0; j < totTask; j++)
			iTasksAll[j] = 0;	

		for (int i = 0; i < gridSize; i++)
		{
			//for (int j = 0; j < totTask; j++)
			{
				//int thisTask = tmpTaskOnGridNode[j][i];
				int thisTask = tmpTaskOnGridNode[i];
	
				if (thisTask > 0)
				{
					//GlobalGrid.globalTaskOnGridNode[i].push_back(thisTask);				
					//iTasksAll[j]++;

					//if (iTasksAll[j] < 10)
					if (thisTask > totTask)
						cout << "node=" << i << " N = " << thisTask << endl;
						//cout << j << ",node=" << i << " N = " << iTasksAll[j] << " " << thisTask << endl;
				}
			}	
		}

#ifdef TEST_BLOCK
		for (int i = 0; i < gridSize; i++)
		{
			GlobalGrid.globalTaskOnGridNode[i].shrink_to_fit();
			int nHalosNode = GlobalGrid.globalTaskOnGridNode[i].size();
	
			if (nHalosNode > 1)		
				cout << "On node = " << i << " N tasks = " << nHalosNode << endl;
		}		
#endif
	}
#endif
};


void Communication::BufferSendRecv()
{
	/* Technical note:
	 * Halo buffers can be allocated directly as vectors.
	 * However, since Particles buffers are allocated as vectors of vectors, the contiguity of ALL the memory blocks is not
	 * guaranteed - thus, it is safer to MPI_Pack all the particle objects into a generic void* buffer and then unpack it
	 */
	vector <Halo> locHalosBufferSend, locHalosBufferRecv;
	void * locPartsBufferSend = NULL; void * locPartsBufferRecv = NULL;
	vector <vector <Particle>> tmpPartsBuffer; 	// Particles will be unpacked here

	size_t locHalosBufferSendSize = 0, locHalosBufferRecvSize = 0;
	size_t locPartsBufferSendSize = 0, locPartsBufferRecvSize = 0;

	int nHalosBufferSend = 1 + locTask, nHalosBufferRecv = 0;
	int bufferPosPartSend=0, bufferPosPartRecv=0;
	int sendTask, recvTask, barrMpi = 0;

	/* Allocate halo buffer */
	locHalosBufferSendSize = nHalosBufferSend * sizeHalo;
	locHalosBufferSend.resize(nHalosBufferSend);

	/* Now compute the size of the particle buffer */
	locPartsBufferSendSize = 0;

	for (int iP = 0; iP < nHalosBufferSend; iP ++)
	{
		locHalosBufferSend[iP] = locHalos[iP];
		locPartsBufferSendSize += locHalos[iP].nPart * sizePart;
	}

	/* Local particles to be sent will be mpi-packed into this buffer */
	locPartsBufferSend = (void *) malloc(locPartsBufferSendSize);

	/* Pack all the selected particles into a single buffer */
	for (int iP = 0; iP < nHalosBufferSend; iP ++)
		MPI_Pack(&locParts[iP][0], locHalos[iP].nPart * sizePart, MPI_BYTE, 
 			  locPartsBufferSend, locPartsBufferSendSize, &bufferPosPartSend, MPI_COMM_WORLD);

	// FIXME this is just a temporary setting
	/* RecvTask is recieving from LocTask and SendTask is sending to LocTask */
	recvTask = (locTask + totTask - 1) % (totTask);
	sendTask = (locTask + 1) % (totTask);

	/* Communicate halo and particle numbers to be sent across tasks */
	MPI_Sendrecv(&nHalosBufferSend, 1, MPI_INT, recvTask, 0, &nHalosBufferRecv, 1, MPI_INT, sendTask, 0, MPI_COMM_WORLD, &status);
	
	/* You have to tell the size of the message */
	MPI_Sendrecv(&locHalosBufferSendSize, sizeof(size_t), MPI_BYTE, recvTask, 0, 
		     &locHalosBufferRecvSize, sizeof(size_t), MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &status);
	
	/* Resize the receiving buffer accordingly */
	//locHalosBufferRecv.resize(locHalosBufferRecvSize);
	locHalosBufferRecv.resize(nHalosBufferRecv);

	barrMpi = MPI_Barrier(MPI_COMM_WORLD);
	
	/* Now get the actual message, the halos */
	MPI_Sendrecv(&locHalosBufferSend[0], locHalosBufferSendSize, MPI_BYTE, recvTask, 0, 
		     &locHalosBufferRecv[0], locHalosBufferRecvSize, MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &status);

	/* Halos have been sent across the tasks, so the send buffer can be cleaned now */
	locHalosBufferSend.clear();
	locHalosBufferSend.shrink_to_fit();

	/* Particle swap across tasks - first communicate the size of the buffer to be delivered */
	MPI_Sendrecv(&locPartsBufferSendSize, sizeof(size_t), MPI_BYTE, recvTask, 0, 
		     &locPartsBufferRecvSize, sizeof(size_t), MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &status);

	//cout << locTask << ") expected: " << locPartsBufferRecvSize << "bytes, sending " << locPartsBufferSendSize << "bytes" << endl;
	
	locPartsBufferRecv = malloc(locPartsBufferRecvSize);

	/* Swap the buffer across the tasks */
	MPI_Sendrecv(locPartsBufferSend, locPartsBufferSendSize, MPI_BYTE, recvTask, 0, 
		     locPartsBufferRecv, locPartsBufferRecvSize, MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &status);

	/* Particles have been sent, clear the buffer */
	free(locPartsBufferSend);

	/* This is a vector <vector <Particles>> - type object and needs to be allocated*/
	tmpPartsBuffer.resize(nHalosBufferRecv);

	/* Unpack the */
	for (int iP = 0; iP < nHalosBufferRecv; iP ++)
	{
		int tmpNParts = locHalosBufferRecv[iP].nPart;
		tmpPartsBuffer[iP].resize(tmpNParts);

		MPI_Unpack(locPartsBufferRecv, locPartsBufferRecvSize, &bufferPosPartRecv, &tmpPartsBuffer[iP][0], 
				tmpNParts * sizePart, MPI_BYTE, MPI_COMM_WORLD);
	}

	if (locTask == 0)
		cout << "Halo and particle Sendrecv has finished, cleaning up the buffers..." << endl;

	/* The buffer has been unpacked into the tmpPartsBuffer vector of vectors, free some space */
	free(locPartsBufferRecv);

	// TODO There will be a loop, all the particle buffers will be appended to some general buffer, the tmpBuffer can be 
	// freed afterwards...
#ifdef TEST_BLOCK	// FIXME clean up the memory...
	for (int iP = 0; iP < nHalosBufferSend; iP ++)	
	{
		tmpPartsBuffer[iP].clear();
		tmpPartsBuffer[iP].shrink_to_fit();
	}

	tmpPartsBuffer.clear();
	tmpPartsBuffer.shrink_to_fit();
	
#endif		// Test Block
};



