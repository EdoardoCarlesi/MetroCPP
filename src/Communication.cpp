/* MPI communication routines that take care of broadcasting, swapping, packing and send/recv of messages */

#include <mpi.h>
#include <iostream>
#include <vector>
#include <algorithm>

#include "utils.h"
#include "global_vars.h"
#include "Grid.h"
#include "Halo.h"
#include "MergerTree.h"
#include "Communication.h"

using namespace std;


/* TODO
void Comm::Optimize()	// Write some function that optimizes the memory distribution among the tasks
{};
*/

#ifndef ZOOM
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
	int nNonZero = 0, thisTask = 0, thisNode = 0;
	
	if (locTask == 0)
		cout << "Broadcasting node information to all tasks for Grid[" << iUseCat << "]." << endl;

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
	int sendTask = 0, recvTask = 0;
	vector<int> allIndex;
	
	if (locTask == 0)
		cout << "Gathering buffer information... " << endl;

	/* Make sure the buffers are clean */
	if (buffIndexNodeHalo.size() > 0)
		buffIndexNodeHalo.clear();

	if (buffIndexSendHalo.size() > 0)
		buffIndexSendHalo.clear();

	/* Resize the buffers - hold a list of halos/nodes to be received from each task */
	buffIndexSendHalo.resize(totTask);
	buffIndexNodeHalo.resize(totTask);

	for (int iT = 0; iT < totTask-1; iT++ )
	{
		sendTask = sendTasks[iT];
		recvTask = recvTasks[iT];
		//sendTask = recvTasks[iT];
		//recvTask = sendTasks[iT];

		//if (iT != locTask)
		{
			sizeSendNode = GlobalGrid[1].buffNodes[sendTask].size();
	
			//cout << " OnTask= " << locTask << ", nNodes: " << sizeSendNode << endl;

			MPI_Sendrecv(&sizeSendNode, 1, MPI_INT, sendTask, 0, 
			  	     &sizeRecvNode, 1, MPI_INT, recvTask, 0, MPI_COMM_WORLD, &status);
			//MPI_Isend(&sizeSendNode, 1, MPI_INT, iT, 0, MPI_COMM_WORLD, &request_send); 
			//MPI_Irecv(&sizeRecvNode, 1, MPI_INT, iT, 0, MPI_COMM_WORLD, &request_recv);
		
			//MPI_Barrier(MPI_COMM_WORLD);

			/* If a task is sending nodes to another task, then it must also receive some nodes 
			 * (not the same number in general) */
			//if (sizeRecvNode != 0)
			{
				//if (sizeRecvNode > 0)
				buffIndexNodeHalo[recvTask].resize(sizeRecvNode);
				MPI_Sendrecv(&GlobalGrid[1].buffNodes[sendTask][0], sizeSendNode, MPI_INT, sendTask, 0, 
			                 &buffIndexNodeHalo[recvTask][0], sizeRecvNode, MPI_INT, recvTask, 0, MPI_COMM_WORLD, &status);

				//if (sizeSendNode > 0)
				//MPI_Isend(&GlobalGrid[1].buffNodes[sendTask][0], sizeSendNode, MPI_INT, sendTask, 0, MPI_COMM_WORLD, &request_send);

				//if (sizeRecvNode > 0)
				//MPI_Irecv(&buffIndexNodeHalo[recvTask][0], sizeRecvNode, MPI_INT, recvTask, 0, MPI_COMM_WORLD, &request_recv);

				//MPI_Send(&GlobalGrid[1].buffNodes[iT][0], sizeSendNode, MPI_INT, iT, 0, MPI_COMM_WORLD); //, &status); 
				//MPI_Recv(&buffIndexNodeHalo[iT][0], sizeRecvNode, MPI_INT, iT, 0, MPI_COMM_WORLD, &status);
				//MPI_Barrier(MPI_COMM_WORLD);
#ifdef VERBOSE
				cout << "Task=" << locTask << " is sending nodes list to=" << recvTask << ", send=" << 
				 	sizeSendNode << ", recv= " << sizeRecvNode << endl;
#endif
				/* Now each task knows which nodes need to be sent and to which task.
				 * We collect the halo indexes corresponding to all of these nodes */
				for (int iN = 0; iN < buffIndexNodeHalo[recvTask].size(); iN++)
				{	
					thisNode = buffIndexNodeHalo[recvTask][iN];
					allIndex = GlobalGrid[1].haloOnGridNode[thisNode];
					
					for (int iH = 0; iH < allIndex.size(); iH++)
					{
						buffIndexSendHalo[recvTask].push_back(allIndex[iH]);

						if (buffIndexSendHalo[recvTask][iH] > locHalos[iUseCat].size())
							cout << "WARNING. On Task=" << locTask << ", toTask=" << recvTask
								<< ", BuffSize=" << buffIndexSendHalo[iT].size() 
								<< ", LocHSize=" << locHalos[iUseCat].size() 
								<< ", IndexSend=" << buffIndexSendHalo[iT][iH] << endl;	
					}
				}

			}
		}
	} 

#ifdef VERBOSE
	cout << "[ TASK = " << locTask << "] has " << buffIndexNodeHalo[recvTask].size() << " on the index list. " << endl;
#endif

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

	/*
	for (int iT = 0; iT < totTask-1; iT++)
	{
			
			int iTsend = sendTasks[iT];  
			int iTrecv = recvTasks[iT]; 
		
			if (iTsend == iTrecv)
				cout << ".... (" << iT << ") WARNING. on task " << locTask << " iSend == iRecv. [" << iTsend<< "] Possible deadlocks." << endl; 
				
	}
	*/

#ifdef VERBOSE
	for (int iT = 0; iT < totTask-1; iT++)
		cout << iT << ") on task= " << locTask << ", send=" << sendTasks[iT] << ", recv=" << recvTasks[iT] << endl; 
#endif
};

#endif	// ifndef ZOOM


#ifdef ZOOM
/* In ZOOM mode we distribute both halo catalogs from task 0 across all other tasks 
 * It is assumed that halo catalogs in zoom simulations are small, and they can be stored on each task at the same time. 
 */
void Communication::BufferSendRecv(void)
{
	void *buffSendParts = nullptr; 
	size_t buffSendSizeParts = 0; 

	if (locTask == 0)
		cout << "Broadcasting all halos from task 0..." << endl;

	/* Here the master task decides how to distribute the halos across the other tasks */
	MPI_Bcast(&nLocHalos[iUseCat], 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (locTask != 0)
	{
		locHalos[iUseCat].resize(nLocHalos[iUseCat]);
		locParts[iUseCat].resize(nLocHalos[iUseCat]);
	}

	MPI_Bcast(&locHalos[iUseCat][0], nLocHalos[iUseCat] * sizeHalo, MPI_BYTE, 0, MPI_COMM_WORLD);

#ifdef VERBOSE
	// Sanity check
	if (locTask != 0)
	{
		cout << locTask << ") has " << locHalos[iUseCat].size() << " halos " << endl;
		locHalos[iUseCat][0].Info();
	}
#endif

	/* Initialize the id2Index map on each task, we are in zoom mode here */
	for (int iH = 0; iH < locHalos[iUseCat].size(); iH++)
	{
		id2Index[locHalos[iUseCat][iH].ID] = iH;
		//if (iH < 5 && locTask == 0)
		//	cout << "====>" << locHalos[iUseCat][iH].ID << " " << id2Index[locHalos[iUseCat][iH].ID] << endl;
	}	

	if (locTask == 0)
		cout << "Done." << endl;
	
	/*
	 * 		COMMUNICATE PARTICLE BUFFERS: Only in full operational mode.
	 * 		When resuming the MTrees we do not need to read/distribute the particles files.
	 */
	if (runMode == 0 || runMode == 2)
	{
#ifdef VERBOSE
		if (locTask == 0)
			cout << "MPI_Packing particle packet to be broadcasted... " << endl;
#endif


		if (locTask == 0)
		{
			int posSendPart = 0;

			for (int iH = 0; iH < nLocHalos[iUseCat]; iH ++)
				for (int iT = 0; iT < nPTypes; iT ++)
					buffSendSizeParts += locHalos[iUseCat][iH].nPart[iT] * sizePart;

			buffSendParts = (void *) malloc(buffSendSizeParts);
		
			/* Pack all the selected particles into a single buffer */
			for (int iH = 0; iH < nLocHalos[iUseCat]; iH ++)
			{
				for (int iP = 0; iP < nPTypes; iP ++)
				{
					if (locHalos[iUseCat][iH].nPart[iP] > 0)
					{
					//	cout << iH << " " << locParts[iUseCat][iH][iP].size() << " " << 
					//			locHalos[iUseCat][iH].nPart[iP] << " " << endl;

						MPI_Pack(&locParts[iUseCat][iH][iP][0], locHalos[iUseCat][iH].nPart[iP] * sizePart, MPI_BYTE, 
							  buffSendParts, buffSendSizeParts, &posSendPart, MPI_COMM_WORLD);

					}
				}
			}

		} // locTask == 0

		/* Communicate the total number of particles */
		MPI_Bcast(&buffSendSizeParts, sizeof(size_t), MPI_BYTE, 0, MPI_COMM_WORLD);

		if (locTask != 0)
			buffSendParts = (void *) malloc(buffSendSizeParts);

		MPI_Bcast(buffSendParts, buffSendSizeParts, MPI_BYTE, 0, MPI_COMM_WORLD);

		if (locTask == 0)
			cout << "Particle buffers broadcasted, unpacking..." << endl;

		if (locTask != 0)
		{
			int posSendPart = 0; int nTmpPart = 0;

			/* Unpack the particle buffer */
			for (int iH = 0; iH < nLocHalos[iUseCat]; iH ++)
			{
				locParts[iUseCat][iH].resize(nPTypes);

				for (int iT = 0; iT < nPTypes; iT ++)
				{	
					nTmpPart = locHalos[iUseCat][iH].nPart[iT];

					if (nTmpPart > 0)
					{
						locParts[iUseCat][iH][iT].resize(nTmpPart);
						MPI_Unpack(buffSendParts, buffSendSizeParts, &posSendPart, &locParts[iUseCat][iH][iT][0], 
								nTmpPart * sizePart, MPI_BYTE, MPI_COMM_WORLD);
					}
				}
			}
		}	// Unpack on tasks != 0

	/* Particles have been unpacked, free the buffer */
	free(buffSendParts);

	if (locTask == 0)
		cout << "Freed buffer, halo and particle have been broadcasted to all tasks. " << endl;

	} // if runMode == 0 or == 2
};


#else	/* ZOOM MODE - using full box halos */


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
	void *buffSendParts = nullptr, *buffRecvParts = nullptr;
	size_t buffSendSizeParts = 0, buffRecvSizeParts = 0;
	size_t buffSendSizeHalos = 0, buffRecvSizeHalos = 0;
	int nTmpPart = 0;

	/* Determine the order of sending and receiving tasks to avoid gridlocks and make it consistent through
	 * all the tasks  */
	SetSendRecvTasks();

	/* First communicate the list of nodes to be sent and received by every task 
	 * The list of nodes also contains the list of haloes associated to them  */
	ExchangeBuffers();

	//MPI_Barrier(MPI_COMM_WORLD);

	if (locTask == 0)
		cout << "Sending and receiving halos in the buffer regions..." << endl; 

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

		if (nBuffSendHalos == 0)
			cout << "Task=" << locTask << " has zero send buffer to " << sendTask << endl;  

		for (int iP = 0; iP < nBuffSendHalos; iP ++)
		{
			int iH = buffIndexSendHalo[sendTask][iP];
			
			buffSendHalos.push_back(locHalos[iUseCat][iH]);
			
			for (int iT = 0; iT < nPTypes; iT ++)
				buffSendSizeParts += locHalos[iUseCat][iH].nPart[iT] * sizePart;
		
			if (iH > locHalos[iUseCat].size())	// Sanity check 
				cout << "WARNING. Halo index " << iH << " not found locally (locHalos). " 
					<< "Required in the send buffer from task=" << locTask << " to task=" << sendTask << endl;
		}

		buffSendSizeHalos = nBuffSendHalos * sizeHalo;

#ifdef VERBOSE
		cout << iT << ") On task=" << locTask << ") sending " << nBuffSendHalos << " halos to " << sendTask <<endl;
		cout << iT << ") On task=" << locTask << ") sending " << buffSendSizeHalos << " b to " << sendTask <<endl;
#endif

		//if (sendTask != recvTask) {
		MPI_Sendrecv(&nBuffSendHalos, 1, MPI_INT, sendTask, 0, 
			     &nBuffRecvHalos, 1, MPI_INT, recvTask, 0, MPI_COMM_WORLD, &status);
		//} else {
		//MPI_Irecv(&nBuffRecvHalos, 1, MPI_INT, recvTask, 0, MPI_COMM_WORLD, &request_recv);
		//MPI_Isend(&nBuffSendHalos, 1, MPI_INT, sendTask, 0, MPI_COMM_WORLD, &request_send);
		//}

		//if (nBuffSendHalos == 0)
		{
			Halo dummyHalo;
		//	buffSendHalos.push_back(dummyHalo); //resize(1);
		//	buffSendSizeHalos = sizeHalo;
		}

		//if (nBuffRecvHalos > 0)
		{
			buffRecvSizeHalos = nBuffRecvHalos * sizeHalo;
			buffRecvHalos.resize(nBuffRecvHalos);	
		//} else {
			//cout << "On task " << locTask << " buffer recv halos is " << nBuffRecvHalos << endl;
		//	Halo dummyHalo;
		//	buffRecvHalos.push_back(dummyHalo); //resize(1);
		//	buffRecvSizeHalos = sizeHalo;
		}

#ifdef VERBOSE
		cout << iT << ") On task=" << locTask << ") recving " << nBuffRecvHalos << " halos from " << recvTask <<endl;
		cout << iT << ") On task=" << locTask << ") recving " << buffRecvSizeHalos << " b from " << recvTask <<endl;
#endif

		//MPI_Barrier(MPI_COMM_WORLD);

	//	if (sendTask != recvTask) {
		MPI_Sendrecv(&buffSendHalos[0], buffSendSizeHalos, MPI_BYTE, sendTask, 0, 
			     &buffRecvHalos[0], buffRecvSizeHalos, MPI_BYTE, recvTask, 0, MPI_COMM_WORLD, &status);

	//	} else {
		/* Non-blocking communication, avoids deadlocks */
		//if (nBuffSendHalos > 0)
		//	MPI_Isend(&buffSendHalos[0], buffSendSizeHalos, MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &request_send);
		//if (nBuffRecvHalos > 0)
		//	MPI_Irecv(&buffRecvHalos[0], buffRecvSizeHalos, MPI_BYTE, recvTask, 0, MPI_COMM_WORLD, &request_recv);

		if (nBuffRecvHalos > 0)
		{
			cout << "OnTask " << locTask << " recv: " << recvTask << " size: " << buffRecvSizeHalos << endl;

			//if (buffRecvSizeHalos > 0) 
		//		MPI_Irecv(&buffRecvHalos[0], buffRecvSizeHalos, MPI_BYTE, recvTask, 0, MPI_COMM_WORLD, &request_recv);
		}
		//}
	
		//}
		/*
		try {
		} catch (MPI::Exception e) {
			cout << "MPI ERROR EXCEPTION ON Task: " << locTask << ", send: " << sendTask 
				<< " recv: " << recvTask << ", sizeS: " << buffSendSizeHalos << ", sizeR: " << buffRecvSizeHalos << endl;
		}
		*/

		/* Add the halos to the local buffer AND a local grid node which is now part of the buffer grid */
		for (int iH = 0; iH < nBuffRecvHalos; iH++)
			locBuffHalos.push_back(buffRecvHalos[iH]);

		/* Clean the recv halo buffer */
		if (buffRecvHalos.size() > 0)
		{
			buffRecvHalos.clear();
			buffRecvHalos.shrink_to_fit();
		}

		/* Clean the send buffer after Sendrecv */
		if (buffSendHalos.size() > 0)
		{
			buffSendHalos.clear();
			buffSendHalos.shrink_to_fit();
		}

		cout << "Halos sent/recv on " << locTask << " in mode " << runMode << endl;

	/*
	 * 		COMMUNICATE PARTICLE BUFFERS
	 */
	
	if (runMode == 0 || runMode == 2)	// Exchange particles only if running in mode 0 or 2	(i.e. building the raw trees from scratch)
	{

		/* Local particles to be sent will be mpi-packed into this buffer */
		buffSendParts = (void *) malloc(buffSendSizeParts);

		locBuffParts.resize(locBuffHalos.size());

#ifdef VERBOSE
		cout << "OnTask=" << locTask <<  ", allocating " 
			<< buffSendSizeParts/1024/1024 << "MB particle send buffer. " << endl;
#else
		if (locTask == 0 && iT == 0)
			cout << "Allocating " << buffSendSizeParts/1024/1024 << " MB for the particle send buffer. " << endl;
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

		//MPI_Isend(&buffSendSizeParts, sizeof(size_t), MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &request_send);
		//MPI_Irecv(&buffRecvSizeParts, sizeof(size_t), MPI_BYTE, recvTask, 0, MPI_COMM_WORLD, &request_recv);

		/* Check that the receiving buffer is a null pointer before allocating it */
		if (buffRecvParts != nullptr)
		{
			free(buffRecvParts);
			buffRecvParts = nullptr;
		}

		buffRecvParts = (void *) malloc(buffRecvSizeParts);

#ifdef VERBOSE
		cout << "Allocating " << buffRecvSizeParts/1024/1024 << "MB recv buffer for " << nBuffRecvHalos << endl;
#else
		if (locTask == 0 && iT == 0)
			cout << "Allocating " << buffRecvSizeParts/1024/1024 
				<< "MB for the particle recv buffer. " << endl;
#endif

		/* Send and recv the MPI_Pack-ed buffers */
		MPI_Sendrecv(buffSendParts, buffSendSizeParts, MPI_BYTE, sendTask, 0, 
			     buffRecvParts, buffRecvSizeParts, MPI_BYTE, recvTask, 0, MPI_COMM_WORLD, &status);

		//if (nBuffSendHalos > 0)
		//	MPI_Isend(buffSendParts, buffSendSizeParts, MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &request_send);
		//if (nBuffRecvHalos > 0)
		//	MPI_Irecv(buffRecvParts, buffRecvSizeParts, MPI_BYTE, recvTask, 0, MPI_COMM_WORLD, &request_recv);

		/* Particles have been sent, free the buffer */
		free(buffSendParts);

		int posRecvPart = 0;

		//cout << "UNPACKING " << nBuffRecvHalos << " on Task " << locTask << endl;  

		/* Unpack the particle buffer */
		for (int iH = 0; iH < nBuffRecvHalos; iH++)
		{
			locBuffParts[iBuffTotHalo].resize(nPTypes);

			for (int iT = 0; iT < nPTypes; iT++)
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
			
			iBuffTotHalo++;
		}
	}	/* runMode == 0 */

}	/* Loop on all the send/recv tasks */

	/* Now assign the halos on the buffer to the respective nodes */
	for (int iH = 0; iH < locBuffHalos.size(); iH++)
		GlobalGrid[iUseCat].AssignToGrid(locBuffHalos[iH].X, -iH-1);	// iH is negative - this is used for halos on the buffer, the -1 is added to avoid overlap with halo num 0

	//MPI_Barrier(MPI_COMM_WORLD);

#ifdef VERBOSE
	cout << "Gathered " << locBuffHalos.size() << " halos in the buffer on task=" << locTask << endl;
	cout << "Gathered " << iBuffTotPart << " parts in the buffer on task=" << locTask << endl;
#else
	if (locTask == 0)
		cout << "Gathered halo and particle buffers. " << endl;
#endif
};

#endif	// ZOOM mode


#ifdef ZOOM	
/* Once the forward correlations of the trees have been built, we communicate orphan halo properties across different tasks */
void Communication::SyncOrphanHalos()
{
	int locOrphans = orphanHaloIndex.size();
	int totOrphans = 0;
	int *posOrphans, *indexOrphans, *dispOrphans;

	if (locTask == 0)
		cout << "Synchronizing orphan halo information across all tasks..." << endl;

	//cout << "-->Task=" << locTask << " has " << locOrphans << endl;

	posOrphans = (int *) calloc(totTask, sizeof(int));
	dispOrphans = (int *) calloc(totTask, sizeof(int));
	
	/* Get the total number of orphan halos across all tasks */
	MPI_Allreduce(&locOrphans, &totOrphans, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	indexOrphans = (int *) calloc(totOrphans, sizeof(int));

	if (locTask == 0)
		cout << "In total, " << totOrphans << " orphan halos have been found." << endl;

	//for (int iL = 0; iL < locOrphans; iL++)
	//	cout << "-->Task=" << locTask << " " << orphanHaloIndex[iL] << " " << iL << endl;

	/* Let every task know how many orphan halo ids it needs to receive from each task */
	MPI_Allgather(&locOrphans, 1, MPI_INT, posOrphans, 1, MPI_INT, MPI_COMM_WORLD);

	//for (int iT = 0; iT < totTask; iT++)
	//	cout << "@ Task=" << locTask << " " << posOrphans[iT] << " " << endl;

	MPI_Barrier(MPI_COMM_WORLD);

	dispOrphans[0] = 0; 

	for (int iD = 0; iD < totTask-1; iD++)
		dispOrphans[iD+1] = posOrphans[iD] + dispOrphans[iD]; 

	/* Let every task know which indexes to receive */
	MPI_Allgatherv(&orphanHaloIndex[0], orphanHaloIndex.size(), 
		MPI_INT, indexOrphans, posOrphans, dispOrphans, MPI_INT, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	//for (int iL = 0; iL < locOrphans; iL++)
	//	cout << "### Task= " << locTask << ", " << iL << ", " << indexOrphans[iL] << endl;
	//for (int iL = 0; iL < totOrphans; iL++)
	//	cout << "### Task= " << locTask << ", " << iL << ", " << indexOrphans[iL] << endl;

	/* Now each task updates the "1" locHalo vector with orphan halos to keep track of them at the next step */
	for (int iL = 0; iL < totOrphans; iL++)
	{
		int nOrphanSteps = 0;
		int thisIndex = indexOrphans[iL];
		unsigned long long int thisID = locHalos[0][thisIndex].ID;
	
		if (facOrphanSteps < 1)
		{
			if (locTask == 0)
				cout << "ERROR. facOrphanSteps not set or set to zero. Check the .cfg file." << endl;  

			exit(0);
		}

		nOrphanSteps = int (locHalos[0][thisIndex].nPart[nPTypes] / facOrphanSteps) + 1;	

		/* Only keep track of the orphans for a number of steps smaller than maxOrphanSteps */
		if (locHalos[0][thisIndex].nOrphanSteps < nOrphanSteps)
		{

			/* When reading the tree files, we do not use particles and locMTrees */
			if (runMode == 0 || runMode == 2)
			{
				locMTrees[0][thisIndex].isOrphan = true;
				locMTrees[0][thisIndex].idProgenitor.push_back(thisID);
				locMTrees[0][thisIndex].indexProgenitor.push_back(nLocHalos[1]);

				locParts[1].push_back(locParts[0][thisIndex]);

				/* Remember to loop also over this halo at the next step 
				 * locTreeIndex keeps track of all the halos to be used on each task in the backward comparison */
				locTreeIndex.push_back(nLocHalos[1]);
			}

			/* The orphan (token) halo is stored in memory for the next step */
			locHalos[1].push_back(locHalos[0][thisIndex]);
			locHalos[1][nLocHalos[1]].nOrphanSteps += 1;
			locHalos[1][nLocHalos[1]].isToken = true;
	
			if (locTask == 0)
				if (locHalos[1][nLocHalos[1]].nOrphanSteps > 1 && nOrphanSteps > 4)
				cout << "Orph=" << locTask << " " << nLocHalos[1] << " " << locHalos[1][nLocHalos[1]].nOrphanSteps 
				<< "/" << nOrphanSteps << " " << locHalos[1][nLocHalos[1]].ID 
				<< " " << locHalos[1][nLocHalos[1]].nPart[1] << endl; 

			if (runMode == 1)
				id2Index[thisID] = nLocHalos[1];

			nLocHalos[1]++;
		}
	}	// if < nOrphanSteps

	// SANITY CHECK
	//cout << "Task=" << locTask << " has " << nLocHalos[1] << " (" << locHalos[1].size() << ") halos including orphans. " << endl; 

	/* Once the orphans have been found, clear this vector */
	orphanHaloIndex.clear();
	orphanHaloIndex.shrink_to_fit();
};

#else	// We use a simpler function to synchronize the indexes in non zoom mode


void Communication::SyncIndex()
{
	/* Clean the map just in case */
	id2Index.clear();
	
	for (int iH = 0; iH < locHalos[iUseCat].size(); iH++)
		id2Index[locHalos[iUseCat][iH].ID] = iH;

	if (iUseCat == 1)
		for (int iH = 0; iH < locBuffHalos.size(); iH++)
			id2Index[locBuffHalos[iH].ID] = -iH-1;	// Add -1 to avoid overlap with index 0

}
#endif		// ZOOM


void Communication::CleanBuffer()
{
	if (buffIndexNodeHalo.size() > 0)
		for (int iN = 0; iN < buffIndexNodeHalo.size(); iN++)
		{
			buffIndexNodeHalo[iN].clear();
			buffIndexNodeHalo[iN].shrink_to_fit();
		}

	if (buffIndexSendHalo.size() > 0)
		for (int iN = 0; iN < buffIndexSendHalo.size(); iN++)
		{
			buffIndexSendHalo[iN].clear();
			buffIndexSendHalo[iN].shrink_to_fit();
		}
}

