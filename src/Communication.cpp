#include <mpi.h>
#include <iostream>
#include <vector>
#include <algorithm>

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
		cout << "Gathering " << gridSize << " nodes to all tasks for grid=" << iUseCat << endl;

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

				//if (thisTask > 0 && iUseCat == 1)
				if (thisTask > 0)
				{
				//	if (iUseCat == 1)
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
	}

	//cout << iNonZero << " duplicate nodes on task=" << locTask << endl;

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
	buffNodeHalo.resize(totTask);
	buffSendHalo.resize(totTask);
	nBuffRecvHalo.resize(totTask);

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
				buffNodeHalo[iT].resize(sizeRecvNode);
				MPI_Sendrecv(&GlobalGrid[1].buffNodes[iT][0], sizeSendNode, MPI_INT, iT, 0, 
				                        &buffNodeHalo[iT][0], sizeRecvNode, MPI_INT, iT, 0, MPI_COMM_WORLD, &status);

				//cout << "Task=" << locTask << " to=" << iT << ", send=" << 
				// sizeSendNode << ", recv= " << sizeRecvNode << endl;

				/* Now each task knows which nodes need to be sent and to which task.
				 * We collect the halo indexes corresponding to all of these nodes */
				for (int iN = 0; iN < buffNodeHalo[iT].size(); iN++)
				{	
					thisNode = buffNodeHalo[iT][iN];
					allIndex = GlobalGrid[1].haloOnGridNode[thisNode];
					
					for (int iH = 0; iH < allIndex.size(); iH++)
					{
						buffSendHalo[iT].push_back(allIndex[iH]);
						if (buffSendHalo[iT][iH] > locHalos[iUseCat].size())
							cout << "WARNING---> " << locTask << ", " << iT << ") " << buffSendHalo[iT][iH] << endl;	
						//cout << iH << ", " << iT << "] " << thisNode << " " << buffSendHalo[iT][iH] << endl;
					}
				}

				sizeSendHalo = buffSendHalo[iT].size();

				/* Communicate how many halos each task will send/recv */
				MPI_Sendrecv(&sizeSendHalo, 1, MPI_INT, iT, 0, 
				        &nBuffRecvHalo[iT], 1, MPI_INT, iT, 0, MPI_COMM_WORLD, &status);

				//if (locTask == 0)	
				cout << "Task=" << locTask << " to=" << iT << " send=" << sizeSendHalo << " halos;"  
					<< " from=" << iT << " recv=" << nBuffRecvHalo[iT] << " halos." << endl;
			}

			//cout << locTask << " from task= " << iT << ", requests n nodes: " << buffNodeHalo[iT].size() << endl; 
		}
	} 

	/* Now every task knows what to send and what to receive from/to every other task */
	if (locTask == 0)
		cout << "Done." << endl;

};


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

	for (int iT = 0; iT < totTask-1; iT++)
		cout << iT << ") on task= " << locTask << ", send=" << sendTasks[iT] << ", recv=" << recvTasks[iT] << endl; 

};


void Communication::BufferSendRecv()
{
	/* 
	 * First communicate the list of nodes to be sent and received by every task 
	 * The list of nodes also contains the list of haloes associated to them, the	
	 */
	ExchangeBuffers();

	/* Determine the order of sending and receiving tasks to avoid gridlocks */
	SetSendRecvTasks();

	MPI_Barrier(MPI_COMM_WORLD);

	/* Technical note:
	 * Halo buffers can be allocated directly as vectors.
	 * However, since Particles buffers are allocated as vectors of vectors, the contiguity of ALL the memory blocks is not
	 * guaranteed - thus, it is safer to MPI_Pack all the particle objects into a generic void* buffer and then unpack it
	 */
	vector <Halo> buffSendHalos, buffRecvHalos;
	vector <vector<vector<unsigned long long int>>> tmpBuffParts; 	

	/* These buffers hold the number of halos to be sent/recvd to/from multiple tasks */
	int buffPosSendPart = 0, buffPosRecvPart = 0;
	int nBuffSendHalos, nBuffRecvHalos, recvTask, sendTask;
	int buffPosSendHalo = 0, buffPosRecvHalo = 0;

	/* Buffers and buffer sizes */
	void *locBuffSendParts = NULL, *locBuffRecvParts = NULL;
	void *locBuffSendHalos = NULL, *locBuffRecvHalos = NULL;
	size_t locBuffSendSizeHalos = 0, locBuffRecvSizeHalos = 0;
	size_t locBuffSendSizeParts = 0, locBuffRecvSizeParts = 0;

//	if (locTask == 0)
//		cout << "Sending and receiving halos across all tasks..." << flush; 

	/* Loop on all the tasks, minus one - the local one */
	//for (int iT = 0; iT < totTask; iT++)
	//for (int iT = 0; iT < 1; iT++)
	for (int iT = 0; iT < totTask-1; iT++)
	{
		sendTask = sendTasks[iT];
		recvTask = recvTasks[iT];
		//sendTask = iT;
		//recvTask = iT;
	
		//if (locTask == 0)
		//if (iT != locTask)
		{
	
			nBuffSendHalos = buffSendHalo[sendTask].size();
			buffSendHalos.resize(nBuffSendHalos);
/*			

			for (int iP = 0; iP < nBuffSendHalos; iP ++)
			{
				int iBuf = buffSendHalo[sendTask][iP];

				if (iBuf > locHalos[iUseCat].size())
					cout << "WARNING---> " << locTask << ", " << iT << ") " << iBuf << endl;	
			}
*/

			//cout << iT << ") task=" << locTask << ") nsend " << buffSendHalo[sendTask].size() << " to " << sendTask 
		 	//	   << " nrecv " << nBuffRecvHalo[recvTask] << " from " << recvTask << endl;

			/* Allocate halo buffers */
			locBuffSendSizeHalos = sizeHalo * nBuffSendHalos;
			locBuffSendHalos = (void *) calloc(sizeHalo, nBuffSendHalos); //locBuffSendSizeHalos); 

			/* Now compute the size of the particle buffer and copy the halos into the buffer to be communicated */
			locBuffSendSizeParts = 0;

			vector<int>::const_iterator max_v_buf = max_element(buffSendHalo[sendTask].begin(), buffSendHalo[sendTask].end());
			cout << "OnTask=" << locTask << " send=" << sendTask << "  max=" << *max_v_buf 
					<< " maxH=" << locHalos[iUseCat].size() << endl; 

			int nBuffDummy = nBuffSendHalos;
			int buf_pos = 0, pack_pos = 0;
			void *dummyBuff; dummyBuff = (void *) malloc(nBuffDummy * sizeHalo);

			for (int iP = 0; iP < nBuffSendHalos; iP ++)
			//for (int iP = 0; iP < nBuffDummy; iP ++)
			{
				int iH = buffSendHalo[sendTask][iP];
				
				buffSendHalos[iP] = locHalos[iUseCat][iH];

				//Halo *dummyHalo; dummyHalo = new Halo();
				//void *dummyData; dummyData = (void *) calloc(1, sizeHalo);
			
				if (iH > locHalos[iUseCat].size())
					cout << "WARNING " << iH << " not in locHalos " << locTask << "-->" << sendTask << endl;

				//MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
				//MPI_Pack(&locHalos[iUseCat][0], sizeHalo, MPI_BYTE, locBuffSendHalos, 
				//	  nBuffDummy * sizeHalo, &buffPosSendHalo, MPI_COMM_WORLD);
				//MPI_Pack(&locHalos[iUseCat][iH], sizeHalo, MPI_BYTE, dummyBuff, 
				//MPI_Pack(&locHalos[iUseCat][0], sizeHalo, MPI_BYTE, dummyBuff, 
				//	  nBuffDummy * sizeHalo, &buf_pos, MPI_COMM_WORLD);
				//MPI_Pack(&dummyHalo->mTot, sizeHalo, MPI_BYTE, locBuffSendHalos, 
				//MPI_Pack(dummyData, sizeHalo, MPI_BYTE, locBuffSendHalos, 
					  //20 * sizeHalo, &buf_pos, MPI_COMM_WORLD);
					  //nBuffDummy * sizeHalo, &buffPosSendHalo, MPI_COMM_WORLD);

				//locBuffSendHalos.push_back(locHalos[iUseCat][iH]);
				//locBuffSendSizeParts += locHalos[iUseCat][iH].nPart[nPTypes] * sizePart;
				//cout << iH << " " << iP << " " << locBuffSendSizeHalos << " " << sizeHalo << endl;
			}

			free(dummyBuff);

#ifdef TEST
			/* Now get the actual message, the halos */
		//	if (locBuffSendSizeHalos > 0) 
		//		MPI_Send(&locBuffSendHalos[0], locBuffSendSizeHalos, MPI_BYTE, sendTask, 0, MPI_COMM_WORLD);
		
		//	if (locBuffRecvSizeHalos > 0) 
		//		MPI_Recv(&locBuffRecvHalos[0], locBuffRecvSizeHalos, MPI_BYTE, recvTask, 0, MPI_COMM_WORLD, &status);

		//	MPI_Sendrecv(&nBuffSendHalos, 1, MPI_INT, sendTask, 0, 
		//		     &nBuffRecvHalos, 1, MPI_INT, recvTask, 0, MPI_COMM_WORLD, &status);

		//	locBuffRecvSizeHalos = sizeHalo * nBuffRecvHalos;
		//	locBuffRecvHalos = malloc(locBuffRecvSizeHalos);

		//	cout << "On task=" << locTask << ") sending " << nBuffSendHalos << " to " << sendTask <<endl;
		//	cout << "On task=" << locTask << ") recving " << nBuffRecvHalos << " by " << recvTask <<endl;
			
		//	int *testSend, *testRecv; 
		//	testSend = new int[nBuffSendHalos];
		//	testRecv = new int[nBuffRecvHalos];

		//	MPI_Sendrecv(&testSend[0], nBuffSendHalos, MPI_INT, sendTask, 0, 
		//		     &testRecv[0], nBuffRecvHalos, MPI_INT, recvTask, 0, MPI_COMM_WORLD, &status);

		//	MPI_Sendrecv(&locBuffSendHalos[0], locBuffSendSizeHalos, MPI_BYTE, recvTask, 0, 
		//		     &locBuffRecvHalos[0], locBuffRecvSizeHalos, MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &status);

		//cout << iT << ") Send " << locBuffSendHalos.size()/sizeHalo << " h on t = " << locTask << "  to  " << sendTask << endl;  
		//cout << iT << ") Recv " << locBuffRecvHalos.size()/sizeHalo << " h on t = " << locTask << " from " << recvTask << endl;  

			/* Halos have been sent across the tasks, so the send buffer can be cleaned now */
			if (locBuffSendSizeHalos > 0) 
			{
				locBuffSendHalos.clear();
				locBuffSendHalos.shrink_to_fit();
			}

			if (locBuffRecvSizeHalos > 0) 
			{
				locBuffRecvHalos.clear();
				locBuffRecvHalos.shrink_to_fit();
			}
#endif
		} /* if locTask != iT */

	}	/* Loop on iT, all the tasks */

	if (locTask == 0)
		cout << "Done." << endl;

#ifdef TEST_BLOCK	
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
	//MPI_Sendrecv(&nHalosBufferSend, 1, MPI_INT, recvTask, 0, &nHalosBufferRecv, 1, MPI_INT, sendTask, 0, MPI_COMM_WORLD, &status);

	//barrMpi = MPI_Barrier(MPI_COMM_WORLD); // Maybe it's useless...

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

	// FIXME clean up the memory...
	tmpPartsBuffer.clear();
	tmpPartsBuffer.shrink_to_fit();
	

#endif		// Test Block
};



