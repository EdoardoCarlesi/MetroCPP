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

/* This function collects all the non-empty grid nodes and broadcasts them to all the tasks.
 * Some nodes may be shared among several tasks, so we need to be careful when MPI_Gather-ing them.
 * In the end, a list of all the nodes of the grid and the tasks that hold them is shared among all the
 * tasks, so that everytime a task needs to access a chunk of the box knows where it needs to look for it. 
 */
void Communication::BroadcastAndGatherGrid()
{
	size_t gridSize = GlobalGrid[iUseCat].nNodes, globalGridSize = 0;
	vector<int> tmpTaskOnGridNode, allNonZeroTasks, allNonZeroNodes;
	int iNonZero = 0, nNonZero = 0, thisTask = 0, thisNode = 0;
	
	if (locTask == 0)
		cout << "Gathering " << gridSize << " nodes to all tasks for grid=" << iUseCat << endl;

	if (locTask == 0)
		tmpTaskOnGridNode.resize(totTask * gridSize);
	
	/* Gather all the nodes on totTask * gridSize array, which can be large, but we need to check for nodes shared between
	 * several processors so first we collect everything separately
	 */
	MPI_Gather(&GlobalGrid[iUseCat].taskOnGridNode[0], gridSize, MPI_INT, 
			&tmpTaskOnGridNode[0], gridSize, MPI_INT, 0, MPI_COMM_WORLD);

	/* Loop on all the huge grid to take care of the nodes which might be assigned to several tasks */
	if (locTask == 0)
	{
		for (int j = 0; j < totTask; j++)
			for (int i = 0; i < gridSize; i++)
			{
				thisTask = tmpTaskOnGridNode[i + gridSize * j];

				//if (thisTask > 0 && iUseCat == 1)
				if (thisTask > 0)
				{
				//	if (iUseCat ==1)
				//		cout << "task=" << thisTask << " " << i << endl;

					allNonZeroTasks.push_back(thisTask);	
					allNonZeroNodes.push_back(i);	
				}
			}	

		nNonZero = allNonZeroTasks.size();
	}	// if (locTask == 0)

	MPI_Bcast(&nNonZero, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Allocate the receiving buffers on all tasks
	if (locTask != 0)
	{
		//cout << "Recving: " << nNonZero << " nodes on task " << locTask << endl;
		allNonZeroTasks.resize(nNonZero);
		allNonZeroNodes.resize(nNonZero);
	}

	MPI_Bcast(&allNonZeroNodes[0], nNonZero, MPI_INT, 0, MPI_COMM_WORLD);	
	MPI_Bcast(&allNonZeroTasks[0], nNonZero, MPI_INT, 0, MPI_COMM_WORLD);	
	
	//if (locTask == 0)
	//cout << "Halo positions on the grid nodes have been broadcasted to all tasks for grid=" << iUseCat << endl;

	/* Now allocate the GLOBAL grid informations and assign the task/node connection on every task */
	GlobalGrid[iUseCat].globalTaskOnGridNode.resize(gridSize);

	for (int i = 0; i < nNonZero; i++)
	{
		thisNode = allNonZeroNodes[i];
		thisTask = allNonZeroTasks[i];
		GlobalGrid[iUseCat].globalTaskOnGridNode[thisNode].push_back(thisTask);

		//cout << i << " " <<  thisNode << " task=" << thisTask << endl;

		// Sanity check // all the nodes are 0 by default
/*		if (GlobalGrid.globalTaskOnGridNode[i].size() > 1)
		{
			int taskOne = GlobalGrid.globalTaskOnGridNode[thisNode][0];
			int taskTwo = GlobalGrid.globalTaskOnGridNode[thisNode][1];
			
			iNonZero++;
			//cout << "Task=" << locTask << " node=" << thisNode << " is shared on " << taskOne << " and " << taskTwo << endl;
		}
*/
	}

	//cout << iNonZero << " duplicate nodes on task=" << locTask << endl;

	// Release some memory from the buffers 
	allNonZeroNodes.clear();	allNonZeroNodes.shrink_to_fit();
	allNonZeroTasks.clear();	allNonZeroTasks.shrink_to_fit();
};


void Communication::BufferSendRecv()
{
	/* Technical note:
	 * Halo buffers can be allocated directly as vectors.
	 * However, since Particles buffers are allocated as vectors of vectors, the contiguity of ALL the memory blocks is not
	 * guaranteed - thus, it is safer to MPI_Pack all the particle objects into a generic void* buffer and then unpack it
	 */
	vector <Halo> locHalosBufferSend, locHalosBufferRecv;
	vector <vector<vector<unsigned long long int>>> tmpPartsBuffer; 	// Particles will be unpacked here

	void * locPartsBufferSend = NULL; void * locPartsBufferRecv = NULL;

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
		locHalosBufferSend[iP] = locHalos[iUseCat][iP];
		locPartsBufferSendSize += locHalos[iUseCat][iP].nPart[nPTypes] * sizePart;
	}

	/* Local particles to be sent will be mpi-packed into this buffer */
	locPartsBufferSend = (void *) malloc(locPartsBufferSendSize);

	/* Pack all the selected particles into a single buffer */
	for (int iP = 0; iP < nHalosBufferSend; iP ++)
		for (int iT = 0; iT < nPTypes; iT ++)
			MPI_Pack(&locParts[iUseCat][iP][iT][0], locHalos[iUseCat][iP].nPart[nPTypes] * sizePart, MPI_BYTE, 
 				  locPartsBufferSend, locPartsBufferSendSize, &bufferPosPartSend, MPI_COMM_WORLD);

	// FIXME this is just a temporary setting: should read from the grid requesting the nodes on the buffer
	/* RecvTask is recieving from LocTask and SendTask is sending to LocTask */
	recvTask = (locTask + totTask - 1) % (totTask);
	sendTask = (locTask + 1) % (totTask);

	/* Communicate halo and particle numbers to be sent across tasks */
	MPI_Sendrecv(&nHalosBufferSend, 1, MPI_INT, recvTask, 0, &nHalosBufferRecv, 1, MPI_INT, sendTask, 0, MPI_COMM_WORLD, &status);
	
	/* You have to tell the size of the message */
	MPI_Sendrecv(&locHalosBufferSendSize, sizeof(size_t), MPI_BYTE, recvTask, 0, 
		     &locHalosBufferRecvSize, sizeof(size_t), MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &status);
	
	/* Resize the receiving buffer accordingly */
	locHalosBufferRecv.resize(nHalosBufferRecv);

	//barrMpi = MPI_Barrier(MPI_COMM_WORLD); // Maybe it's useless...
	
	/* Now get the actual message, the halos */
	MPI_Sendrecv(&locHalosBufferSend[0], locHalosBufferSendSize, MPI_BYTE, recvTask, 0, 
		     &locHalosBufferRecv[0], locHalosBufferRecvSize, MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &status);

	/* Halos have been sent across the tasks, so the send buffer can be cleaned now */
	locHalosBufferSend.clear();
	locHalosBufferSend.shrink_to_fit();

#ifdef TEST_BLOCK	// FIXME clean up the memory...
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

	for (int iP = 0; iP < nHalosBufferRecv; iP ++)
		tmpPartsBuffer[iP].resize(nPTypes);
	
	/* Unpack the particle buffer */
	for (int iP = 0; iP < nHalosBufferRecv; iP ++)
	{
		for (int iT = 0; iT < nHalosBufferRecv; iT ++)
		{
			int tmpNParts = locHalosBufferRecv[iP].nPart[nPTypes];
			//tmpPartsBuffer[iP].resize(tmpNParts);
			MPI_Unpack(locPartsBufferRecv, locPartsBufferRecvSize, &bufferPosPartRecv, &tmpPartsBuffer[iP][0], 
					tmpNParts * sizePart, MPI_BYTE, MPI_COMM_WORLD);
		}
	}

	if (locTask == 0)
		cout << "Halo and particle Sendrecv has finished, cleaning up the buffers..." << endl;

	/* The buffer has been unpacked into the tmpPartsBuffer vector of vectors, free some space */
	free(locPartsBufferRecv);

	// TODO There will be a loop, all the particle buffers will be appended to some general buffer, the tmpBuffer can be 
	// freed afterwards...
	for (int iP = 0; iP < nHalosBufferSend; iP ++)	
	{
		tmpPartsBuffer[iP].clear();
		tmpPartsBuffer[iP].shrink_to_fit();
	}

	tmpPartsBuffer.clear();
	tmpPartsBuffer.shrink_to_fit();
	
#endif		// Test Block
};



