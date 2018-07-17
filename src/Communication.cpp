/* MPI communication routines that take care of broadcasting, swapping, packing and send/recv of messages */

#include <mpi.h>
#include <iostream>
#include <vector>
#include <algorithm>

#include "utils.h"
#include "global_vars.h"
#include "Grid.h"
#include "Halo.h"
#include "Communication.h"

using namespace std;


/* TODO
void Comm::Optimize()	// Write some function that optimizes the memory distribution among the tasks
{};
*/

/* 
 * This function collects all the non-empty grid nodes and broadcasts them to all the tasks.
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
		cout << "Broadcasting node informations to all tasks for Grid[" << iUseCat << "]." << endl;

	if (locTask == 0)
		tmpTaskOnGridNode.resize(totTask * gridSize);
	
	/* 
	 * Gather all the nodes on totTask * gridSize array, which can be large, but we need to check for nodes shared between
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

				if (thisTask > 0)
				{
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
	
	if (locTask == 0)
		cout << "Halo positions on the grid nodes have been broadcasted to all tasks." << endl;

	/* Now allocate the GLOBAL grid informations and assign the task/node connection on every task */
	GlobalGrid[iUseCat].globalTaskOnGridNode.resize(gridSize);

	for (int i = 0; i < nNonZero; i++)
	{
		thisNode = allNonZeroNodes[i];
		thisTask = allNonZeroTasks[i];
		GlobalGrid[iUseCat].globalTaskOnGridNode[thisNode].push_back(thisTask);
	}

#ifdef VERBOSE
	cout << iNonZero << " duplicate nodes on task=" << locTask << endl;
#endif
	// Release some memory from the buffers 
	allNonZeroNodes.clear();	allNonZeroNodes.shrink_to_fit();
	allNonZeroTasks.clear();	allNonZeroTasks.shrink_to_fit();
};


/*
 * This function communicates the list of nodes that each task needs to send (recv) from every other task
 * Every node keeps track of its nearby halos, these halos are then MPI_Packed and delivered to the tasks 
 * that request them.
 */
void Communication::ExchangeBuffers()
{
	int sizeSendNode = 0, sizeRecvNode = 0, thisNode = 0, sizeSendHalo = 0;
	vector<int> allIndex;
	
	if (locTask == 0)
		cout << "Exchanging buffer among tasks... " << endl;

	// Resize the buffers - hold a list of halos/nodes to be received from each task
	buffIndexNodeHalo.resize(totTask);
	buffIndexSendHalo.resize(totTask);

	for (int iT = 0; iT < totTask; iT++ )
	{
		if (iT != locTask)
		{
			sizeSendNode = GlobalGrid[1].buffNodes[iT].size();

			MPI_Sendrecv(&sizeSendNode, 1, MPI_INT, iT, 0, 
			  	     &sizeRecvNode, 1, MPI_INT, iT, 0, MPI_COMM_WORLD, &status);
		
			/* If a task is sending nodes to another task, then it must also receive some nodes 
			 * (not the same number in general) */
			if (sizeRecvNode != 0)
			{
				buffIndexNodeHalo[iT].resize(sizeRecvNode);
				MPI_Sendrecv(&GlobalGrid[1].buffNodes[iT][0], sizeSendNode, MPI_INT, iT, 0, 
				                   &buffIndexNodeHalo[iT][0], sizeRecvNode, MPI_INT, iT, 0, MPI_COMM_WORLD, &status);
#ifdef VERBOSE
				cout << "Task=" << locTask << " to=" << iT << ", send=" << 
				 	sizeSendNode << ", recv= " << sizeRecvNode << endl;
#endif
				/* Now each task knows which nodes need to be sent and to which task.
				 * We collect the halo indexes corresponding to all of these nodes */
				for (int iN = 0; iN < buffIndexNodeHalo[iT].size(); iN++)
				{	
					thisNode = buffIndexNodeHalo[iT][iN];
					allIndex = GlobalGrid[1].haloOnGridNode[thisNode];
					
					for (int iH = 0; iH < allIndex.size(); iH++)
					{
						buffIndexSendHalo[iT].push_back(allIndex[iH]);
						if (buffIndexSendHalo[iT][iH] > locHalos[iUseCat].size())
							cout << "WARNING---> " << locTask << ", " << iT 
								<< ") " << buffIndexSendHalo[iT][iH] << endl;	
					}
				}

			}
		}
	} 

	/* Now every task knows what to send and what to receive from/to every other task */
	if (locTask == 0)
		cout << "Done." << endl;

};


/* Each task is sending and receiving to/from multiple task. 
 * Send/recv tasks need to be defined consistently everywhere to ensure correct communication. */
void Communication::SetSendRecvTasks()
{
	int iTrecv = 0, iTsend = 0;
	int jTrecv = 0, jTsend = 0;

	sendTasks.resize(totTask - 1);
	recvTasks.resize(totTask - 1);

	for (int iT = 0; iT < totTask; iT++)
	{
		iTrecv = (locTask - 1 - iT + totTask) % totTask;
		iTsend = (iT + locTask + 1 + totTask) % totTask;

		if (iTrecv != locTask)
		{
			recvTasks[jTrecv] = iTrecv;
			jTrecv++;
		}		

		if (iTsend != locTask)
		{
			sendTasks[jTsend] = iTsend;
			jTsend++;
		}		
		
	}

#ifdef VERBOSE
	for (int iT = 0; iT < totTask-1; iT++)
		cout << iT << ") on task= " << locTask << ", send=" << sendTasks[iT] << ", recv=" << recvTasks[iT] << endl; 
#endif
};


/* This function first determines the size of the buffers to be broadcasted, then communicates it to all the tasks.
 * Each task packs all the halos and particles that are requested by other halos for comparison into several buffers,
 * which are communicated with a call to MPI_Sendrecv. The buffers are then unpacked into the locBuffHalos and locBuffParts
 * vectors, which contain all of the halos/particles received from the neighbouring subvolumes. */
void Communication::BufferSendRecv()
{
	/* Tasks receiving and sending messages */
	int recvTask = 0, sendTask = 0;
	
	/* Keep track of all the halos & particles in the buffer region */
	int iBuffTotHalo = 0;
	int iBuffTotPart = 0;

	/* These buffers hold the number of halos to be sent/recvd to/from multiple tasks */
	int nBuffSendHalos = 0, nBuffRecvHalos = 0; 
	int buffIndexHalo = 0, buffIndexPart = 0;

	/* Technical note:
	 * Halo buffers can be allocated directly as vectors.
	 * However, since Particles buffers are allocated as vectors of vectors, the contiguity of ALL the memory blocks is not
	 * guaranteed - thus, it is safer to MPI_Pack all the particle objects into a generic void* buffer and then unpack it  */
	vector<Halo> buffSendHalos;
	vector<Halo> buffRecvHalos;

	/* Buffers and buffer sizes */
	void *buffSendParts = NULL, *buffRecvParts = NULL;
	size_t buffSendSizeParts = 0, buffRecvSizeParts = 0;
	size_t buffSendSizeHalos = 0, buffRecvSizeHalos = 0;
	int nTmpPart = 0;

	/* First communicate the list of nodes to be sent and received by every task 
	 * The list of nodes also contains the list of haloes associated to them  */
	ExchangeBuffers();

	/* Determine the order of sending and receiving tasks to avoid gridlocks and make it consistent through
	 * all the tasks  */
	SetSendRecvTasks();

	MPI_Barrier(MPI_COMM_WORLD);

	if (locTask == 0)
		cout << "Exchanging halos in the buffer region across all tasks..." << endl; 

	/*
	 * 		FIRST EXCHANGE HALO BUFFERS
	 */

	/* Loop on all the tasks, minus one - the local one */
	for (int iT = 0; iT < totTask-1; iT++)
	{
		sendTask = sendTasks[iT];
		recvTask = recvTasks[iT];

		/* Now compute the size of the particle buffer and copy the halos into the buffer to be communicated */
		buffSendSizeParts = 0;

		/* At this step, locTask will send nBuffSendHalos to sendTask */
		nBuffSendHalos = buffIndexSendHalo[sendTask].size();

		for (int iP = 0; iP < nBuffSendHalos; iP ++)
		{
			int iH = buffIndexSendHalo[sendTask][iP];
			
			buffSendHalos.push_back(locHalos[iUseCat][iH]);
			
			for (int iT = 0; iT < nPTypes; iT ++)
				buffSendSizeParts += locHalos[iUseCat][iH].nPart[iT] * sizePart;
		
			if (iH > locHalos[iUseCat].size())	// Sanity check 
				cout << "WARNING " << iH << " not in locHalos " << locTask << "-->" << sendTask << endl;
		}

		buffSendSizeHalos = nBuffSendHalos * sizeHalo;

		MPI_Sendrecv(&nBuffSendHalos, 1, MPI_INT, sendTask, 0, 
			     &nBuffRecvHalos, 1, MPI_INT, recvTask, 0, MPI_COMM_WORLD, &status);

		buffRecvSizeHalos = nBuffRecvHalos * sizeHalo;
		buffRecvHalos.resize(nBuffRecvHalos);	

#ifdef VERBOSE
		cout << iT << ") On task=" << locTask << ") sending " << buffSendSizeHalos << " to " << recvTask <<endl;
		cout << iT << ") On task=" << locTask << ") recving " << buffRecvSizeHalos << " by " << sendTask <<endl;
#endif

		MPI_Sendrecv(&buffSendHalos[0], buffSendSizeHalos, MPI_BYTE, sendTask, 0, 
			     &buffRecvHalos[0], buffRecvSizeHalos, MPI_BYTE, recvTask, 0, MPI_COMM_WORLD, &status);

		/* Clean the send buffer after Sendrecv */
		if (buffSendHalos.size() > 0)
		{
			buffSendHalos.clear();
			buffSendHalos.shrink_to_fit();
		}

		/* Add the halos to the local buffer AND a local grid node which is now part of the buffer grid */
		for (int iH = 0; iH < nBuffRecvHalos; iH++)
		{
			locBuffHalos.push_back(buffRecvHalos[iH]);
			buffIndexHalo = locBuffHalos.size() - 1; 
			BufferGrid.AssignToGrid(buffRecvHalos[iH].X, buffIndexHalo);
		}

		/* Clean the recv halo buffer */
		if (buffRecvHalos.size() > 0)
		{
			buffRecvHalos.clear();
			buffRecvHalos.shrink_to_fit();
		}

		/*
		 * 		COMMUNICATE PARTICLE BUFFERS
		 */

		/* Local particles to be sent will be mpi-packed into this buffer */
		buffSendParts = (void *) malloc(buffSendSizeParts);

		locBuffParts.resize(locBuffHalos.size());

#ifdef VERBOSE
		cout << "Allocating " << buffSendSizeParts/1024/1024 << "MB send buffer for " << nBuffSendHalos << endl;
#else
		if (locTask == 0 && iT == 0)
			cout << "Allocating particle send buffer... " << endl;
#endif

		/* MPI_Pack variables - for some reason if defined at the beginning MPI_Pack crashes... */
		int posSendPart = 0;

		/* Pack all the selected particles into a single buffer */
		for (int iI = 0; iI < nBuffSendHalos; iI ++)
		{
			int iH = buffIndexSendHalo[sendTask][iI];

			for (int iP = 0; iP < nPTypes; iP ++)
			{
				if (locHalos[iUseCat][iH].nPart[iP] > 0)
				{
					MPI_Pack(&locParts[iUseCat][iH][iP][0], locHalos[iUseCat][iH].nPart[iP] * sizePart, MPI_BYTE, 
						  buffSendParts, buffSendSizeParts, &posSendPart, MPI_COMM_WORLD);
				}
			}
		}

		/* Communicate the MPI_Pack-ed buffer sizes */
		MPI_Sendrecv(&buffSendSizeParts, sizeof(size_t), MPI_BYTE, sendTask, 0, 
			     &buffRecvSizeParts, sizeof(size_t), MPI_BYTE, recvTask, 0, MPI_COMM_WORLD, &status);

		/* Check that the receiving buffer is NULL before allocating it */
		if (buffRecvParts != NULL)
		{
			free(buffRecvParts);
			buffRecvParts = NULL;
		}

		buffRecvParts = (void *) malloc(buffRecvSizeParts);

#ifdef VERBOSE
		cout << "Allocating " << buffRecvSizeParts/1024/1024 << "MB recv buffer for " << nBuffRecvHalos << endl;
#else
		if (locTask == 0 && iT == 0)
			cout << "Allocating particle recv buffer..." << endl;
#endif

		/* Send and recv the MPI_Pack-ed buffers */
		MPI_Sendrecv(buffSendParts, buffSendSizeParts, MPI_BYTE, sendTask, 0, 
			     buffRecvParts, buffRecvSizeParts, MPI_BYTE, recvTask, 0, MPI_COMM_WORLD, &status);

		/* Particles have been sent, free the buffer */
		free(buffSendParts);

		int posRecvPart = 0;

		/* Unpack the particle buffer */
		for (int iH = 0; iH < nBuffRecvHalos; iH ++)
		{
			locBuffParts[iBuffTotHalo].resize(nPTypes);

			for (int iT = 0; iT < nPTypes; iT ++)
			{	
				nTmpPart = locBuffHalos[iBuffTotHalo].nPart[iT];
				iBuffTotPart += nTmpPart;

				if (nTmpPart > 0)
				{
					locBuffParts[iBuffTotHalo][iT].resize(nTmpPart);
	
					MPI_Unpack(buffRecvParts, buffRecvSizeParts, &posRecvPart, &locBuffParts[iBuffTotHalo][iT][0], 
							nTmpPart * sizePart, MPI_BYTE, MPI_COMM_WORLD);
				}
			
			}

			iBuffTotHalo++;		// Keep track of the total number of halos in the buffer
		}
	}	/* Loop on all the send/recv tasks */

#ifdef VERBOSE
	cout << "Gathered " << locBuffHalos.size() << " halos in the buffer on task=" << locTask << endl;
	cout << "Gathered " << iBuffTotHalo << " halos in the buffer on task=" << locTask << endl;
	cout << "Gathered " << iBuffTotPart << " parts in the buffer on task=" << locTask << endl;
#else
	if (locTask == 0)
		cout << "Gathered halo and particle buffers. " << endl;
#endif
};



