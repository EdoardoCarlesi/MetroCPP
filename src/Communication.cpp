/*
 *   METROC++: MErger TRees On C++, a scalable code for the computation of merger trees in cosmological simulations.
 *   Copyright (C) Edoardo Carlesi 2018-2019
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/* Communication.cpp : 
   MPI communication routines that take care of broadcasting, swapping, packing and send/recv of messages */

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


Communication::~Communication()
{
	sendTasks.clear();
	sendTasks.shrink_to_fit();
	recvTasks.clear();
	recvTasks.shrink_to_fit();

	for (int iN = 0; iN < buffIndexNodeHalo.size(); iN++)
		buffIndexNodeHalo.clear();

	for (int iN = 0; iN < buffIndexSendHalo.size(); iN++)
		buffIndexSendHalo.clear();

	buffIndexSendHalo.clear();
	buffIndexSendHalo.shrink_to_fit();

	buffIndexNodeHalo.clear();
	buffIndexNodeHalo.shrink_to_fit();
};


/* NOTE: In ZOOM mode there is no buffer to be communicated, as everything works on shared memory */

#ifndef ZOOM
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

#ifdef VERBOSE
		if (nBuffSendHalos == 0)
			cout << "Task=" << locTask << " has zero send buffer to " << sendTask << endl;  
#endif

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
		MPI_Sendrecv(&nBuffSendHalos, 1, MPI_INT, sendTask, 0, 
			     &nBuffRecvHalos, 1, MPI_INT, recvTask, 0, MPI_COMM_WORLD, &status);

		buffRecvSizeHalos = nBuffRecvHalos * sizeHalo;
		buffRecvHalos.resize(nBuffRecvHalos);	

#ifdef VERBOSE
		cout << iT << ") On task=" << locTask << ") recving " << nBuffRecvHalos << " halos from " << recvTask <<endl;
		cout << iT << ") On task=" << locTask << ") recving " << buffRecvSizeHalos << " b from " << recvTask <<endl;
#endif
		MPI_Sendrecv(&buffSendHalos[0], buffSendSizeHalos, MPI_BYTE, sendTask, 0, 
			     &buffRecvHalos[0], buffRecvSizeHalos, MPI_BYTE, recvTask, 0, MPI_COMM_WORLD, &status);

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

	/*
	 * 		MPI_Pack & exchange particle buffers
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
			cout << "Allocating " << buffSendSizeParts/1024/1024 << "MB for the particle send buffer. " << endl;
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

		/* Particles have been sent, free the buffer */
		free(buffSendParts);

		int posRecvPart = 0;

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

					for (auto const& partID : locBuffParts[iBuffTotHalo][iT])
					{
				       		Particle thisParticle;
   	                      	        	thisParticle.haloID = locBuffHalos[iBuffTotHalo].ID;
        	                        	thisParticle.type   = iT;
                	                	locMapParts[iUseCat][partID].push_back(thisParticle);
					}
				}
			}
			
			iBuffTotHalo++;
		}
	}	/* runMode == 0 */

}	/* Loop on all the send/recv tasks */

	/* Now assign the halos on the buffer to the respective nodes */
	for (int iH = 0; iH < locBuffHalos.size(); iH++)
		GlobalGrid[iUseCat].AssignToGrid(locBuffHalos[iH].X, -iH-1);	// iH is negative - this is used for halos on the buffer, the -1 is added to avoid overlap with halo num 0

#ifdef VERBOSE
	cout << "Gathered " << locBuffHalos.size() << " halos in the buffer on task=" << locTask << endl;
	cout << "Gathered " << iBuffTotPart << " parts in the buffer on task=" << locTask << endl;
#else
	if (locTask == 0)
		cout << "Gathered halo and particle buffers. " << endl;
#endif
};



/* Once the merger trees are built, we need to communicate the results 
 * 
 */
void Communication::SyncMergerTreeBuffer()
{
	/* Tasks receiving and sending messages */
	int recvTask = 0, sendTask = 0;

	/* These buffers hold the number of halos to be sent/recvd to/from multiple tasks */
	int iBuffTotHaloIDs = 0, nBuffSendHaloIDs = 0, nBuffRecvHaloIDs = 0, buffIndexHaloIDs = 0, nTotBuffHaloIDs = 0;
	
	size_t buffSendSizeHaloIDs = 0, buffRecvSizeHaloIDs = 0;

	/* These vectors contain the main halo and its main progenitor IDs, like [ID1, progID1, ID2, progID2... ]. 
	 * The descendant ID is the i-th index while the progenitor is the i-th+1. */
	vector <uint64_t> buffSendHaloIDs;
	vector <uint64_t> buffRecvHaloIDs;
	vector <uint64_t> totBuffRecvHaloIDs;

	// FIXME maybe these things are already set, no need to call these functions once again
	/* Determine the order of sending and receiving tasks to avoid gridlocks and make it consistent through
	 * all the tasks  */
	SetSendRecvTasks();

	/* First communicate the list of nodes to be sent and received by every task 
	 * The list of nodes also contains the list of haloes associated to them  */
	//ExchangeBuffers();

	if (locTask == 0)
		cout << "Synchronizing merger trees in the buffer regions..." << endl; 

	//for (int iT = 0; iT < totTask-1; iT++)
	//	cout << iT << " on Task " << locTask << " send: " << sendTasks[iT] << " recv: " << recvTasks[iT] << endl;  


	/* Loop on all the tasks except the local one */
	for (int iT = 0; iT < totTask-1; iT++)
	{
		sendTask = sendTasks[iT];
		recvTask = recvTasks[iT];

		/* At this step, locTask will send nBuffSendHalos to sendTask */
		nBuffSendHaloIDs = buffIndexSendHalo[sendTask].size();

#ifdef VERBOSE
		if (nBuffSendHaloIDs == 0)
			cout << "Task=" << locTask << " has zero send buffer to " << sendTask << endl;  
#endif

		buffSendHaloIDs.resize(2 * nBuffSendHaloIDs);

		for (int iP = 0; iP < nBuffSendHaloIDs; iP ++)
		{
			int iH = buffIndexSendHalo[sendTask][iP];
			
			/* i-th and i+1-th IDs in the vector are a descendant/progenitor pair */
			buffSendHaloIDs[2 * iP] = locMTrees[iUseCat][iH].mainHalo.ID;

			/* If the halo on the buffer is orphan, then we assign it its own ID as progenitor */
			if (locMTrees[iUseCat][iH].isOrphan == true)
				buffSendHaloIDs[2 * iP + 1] = locMTrees[1][iH].mainHalo.ID;
			else
				buffSendHaloIDs[2 * iP + 1] = locMTrees[1][iH].idProgenitor[0];
			
			if (iH > locHalos[iUseCat].size())	// Sanity check 
				cout << "WARNING. Halo index " << iH << " not found locally (locHalos). " 
					<< "Required in the send buffer from task=" << locTask << " to task=" << sendTask << endl;
		}

#ifdef VERBOSE
		cout << iT << ") On task=" << locTask << ") sending " << nBuffSendHaloIDs << " halos to " << sendTask <<endl;
		cout << iT << ") On task=" << locTask << ") sending " << buffSendSizeHaloIDs << " b to " << sendTask <<endl;
#endif
		
		/* Communicate the number of halo IDs on the buffer (progenitor only) */
		MPI_Sendrecv(&nBuffSendHaloIDs, 1, MPI_INT, sendTask, 0, 
			     &nBuffRecvHaloIDs, 1, MPI_INT, recvTask, 0, MPI_COMM_WORLD, &status);

		/* Multiply by two to take into account progenitor and descendant */
		buffSendSizeHaloIDs = 2 * nBuffSendHaloIDs * sizeof(uint64_t);
		buffRecvSizeHaloIDs = 2 * nBuffRecvHaloIDs * sizeof(uint64_t);

		nTotBuffHaloIDs += nBuffRecvHaloIDs;

		buffRecvHaloIDs.resize(2 * nBuffRecvHaloIDs);	

#ifdef VERBOSE
		cout << iT << ") On task=" << locTask << ") recving " << nBuffRecvHaloIDs << " halos from " << recvTask <<endl;
		cout << iT << ") On task=" << locTask << ") recving " << buffRecvSizeHaloIDs << " b from " << recvTask <<endl;
#endif
		MPI_Sendrecv(&buffSendHaloIDs[0], buffSendSizeHaloIDs, MPI_BYTE, sendTask, 0, 
			     &buffRecvHaloIDs[0], buffRecvSizeHaloIDs, MPI_BYTE, recvTask, 0, MPI_COMM_WORLD, &status);

		/* Append halo ids to the total recv buffer */
		for (int iH = 0; iH < 2 * nBuffRecvHaloIDs; iH++)
			totBuffRecvHaloIDs.push_back(buffRecvHaloIDs[iH]);

		/* Clean the recv halo buffer */
		if (buffRecvHaloIDs.size() > 0)
		{
			buffRecvHaloIDs.clear();
			buffRecvHaloIDs.shrink_to_fit();
		}

		/* Clean the send buffer after Sendrecv */
		if (buffSendHaloIDs.size() > 0)
		{
			buffSendHaloIDs.clear();
			buffSendHaloIDs.shrink_to_fit();
		}

	}	// Loop on the send/recv tasks

	// TODO fix this map trees in non CMP mode

	/* Now synchronize the connections of the halos on the buffer */
	for (int iH = 0; iH < locBuffHalos.size(); iH++)
	{	
		uint64_t descID = totBuffRecvHaloIDs[2 * iH];
		uint64_t progID = totBuffRecvHaloIDs[2 * iH + 1];

		int thisTreeIndex = thisMapTrees[descID];

		/* Check whether the most likely descendant of the halo in the buffer is the same as the one found in its "original" task. */
		if (locMTrees[1][thisTreeIndex].mainHalo.ID == descID && locMTrees[1][thisTreeIndex].idProgenitor.size() > 0)
		{

			/* If the progenitor with the highest merit is not on this task, replace the (local) highest merit ID with this one.
			 * In this way, this halo now knows that its most likely descendant is not located on this task, and when cleaning the 
			 * connection will not be taken among the progenitors of the old highest merit ID-halo.	*/
			if (locMTrees[1][thisTreeIndex].idProgenitor[0] != progID)
				locMTrees[1][thisTreeIndex].idProgenitor[0] = progID;
		}
	}
}


/* This function collects all the non-empty grid nodes and broadcasts them to all the tasks.
 * Some nodes may be shared among several tasks, so we need to be careful when MPI_Gather-ing them.
 * In the end, a list of all the nodes of the grid and the tasks that hold them is shared among all the
 * tasks, so that everytime a task needs to access a chunk of the box knows where it needs to look for it. */
void Communication::BroadcastAndGatherGrid()
{
	size_t gridSize = GlobalGrid[iUseCat].nNodes, globalGridSize = 0;
	vector<int> tmpTaskOnGridNode, allNonZeroTasks, allNonZeroNodes;
	int nNonZero = 0, thisTask = 0, thisNode = 0;
	
	if (locTask == 0)
		cout << "Exchanging node information among tasks..." << endl;

	if (locTask == 0)
		tmpTaskOnGridNode.resize(totTask * gridSize);
	
	/* Gather all the nodes on totTask * gridSize array, which can be large, but we need to check for nodes shared between
	 * several processors so first we collect everything separately. */
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


/* This function communicates the list of nodes that each task needs to send (recv) from every other task
 * Every node keeps track of its nearby halos, these halos are then MPI_Packed and delivered to the tasks 
 * that request them. */
void Communication::ExchangeBuffers()
{
	int sizeSendNode = 0, sizeRecvNode = 0, thisNode = 0, sizeSendHalo = 0;
	int sendTask = 0, recvTask = 0;
	vector<int> allIndex;
	
	if (locTask == 0)
		cout << "Gathering buffer information..." << flush;

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

		{
			sizeSendNode = GlobalGrid[1].buffNodes[sendTask].size();

			MPI_Sendrecv(&sizeSendNode, 1, MPI_INT, sendTask, 0, 
			  	     &sizeRecvNode, 1, MPI_INT, recvTask, 0, MPI_COMM_WORLD, &status);

			/* If a task is sending nodes to another task, then it must also receive some nodes 
			 * (not the same number in general) */
			{
				buffIndexNodeHalo[recvTask].resize(sizeRecvNode);
				MPI_Sendrecv(&GlobalGrid[1].buffNodes[sendTask][0], sizeSendNode, MPI_INT, sendTask, 0, 
			                 &buffIndexNodeHalo[recvTask][0], sizeRecvNode, MPI_INT, recvTask, 0, MPI_COMM_WORLD, &status);

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
		cout << "done." << endl;

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

#endif		// Non zoom mode (ifdef ZOOM, else, endif)


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

