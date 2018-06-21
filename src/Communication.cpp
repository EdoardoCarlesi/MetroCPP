#include <mpi.h>
#include <iostream>
#include <vector>

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
void Comm::Optimize()	// Write some function that 
{};
*/

void Communication::BufferSendRecv()
{
	int nHalosBufferSend = 1 + locTask, nHalosBufferRecv = 0;
	int bufferPosHaloSend=0, bufferPosHaloRecv=0;
	int bufferPosPartSend=0, bufferPosPartRecv=0;
	int sendTask, recvTask;
	int barrMpi = 0;

	/* Vectors that will store haloes sent from the neighbouring tasks */
	vector<Halo> locHalosBufferNew;
	vector<vector <Particle>> locPartsBufferNew;

	locHalosBufferSendSize = nHalosBufferSend * sizeHalo;
	locHalosBufferSend.resize(nHalosBufferSend);
	locPartsBufferSendSize = 0;

	for (int iP = 0; iP < nHalosBufferSend; iP ++)
	{
		locHalosBufferSend[iP] = locHalos[iP];
		locPartsBufferSendSize += locHalos[iP].nPart * sizePart;
	}

	/*
	try 
	{
		locPartsBufferSend = (void *) malloc(locPartsBufferSendSize);
	} catch (bad_alloc& ba) {
		cout << "Could not allocate " << (locPartsBufferSendSize / 1024 / 1024) << " MB of parts on task=" << locTask << endl; 
	}
	*/

	// FIXME this is just a temporary setting
	/* RecvTask is recieving from LocTask and SendTask is sending to LocTask */
	recvTask = (locTask + totTask - 1) % (totTask);
	sendTask = (locTask + 1) % (totTask);

	cout << "On task=" << locTask << " sending " << locHalosBufferSendSize << " bytes of halos to task=" << recvTask << endl;
	//cout << "On task=" << locTask << " MPI_Pack-ing " << locPartsBufferSend.size() << " bytes of parts to task=" << recvTask << endl;
 
	/* Communicate halo and particle numbers to be sent across tasks */
	MPI_Sendrecv(&nHalosBufferSend, 1, MPI_INT, recvTask, 0, &nHalosBufferRecv, 1, MPI_INT, sendTask, 0, MPI_COMM_WORLD, &status);

	MPI_Sendrecv(&locHalosBufferSendSize, sizeof(size_t), MPI_BYTE, recvTask, 0, 
		     &locHalosBufferRecvSize, sizeof(size_t), MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &status);

#ifdef TEST_BLOCK
	cout << "On task=" << locTask << ", expecting =" << nHalosBufferRecv << " halos from task=" << sendTask << endl;

	MPI_Sendrecv(&locPartsBufferSendSize, sizeof(size_t), MPI_BYTE, recvTask, 0, 
		     &locPartsBufferRecvSize, sizeof(size_t), MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &status);

	
	locHalosBufferRecv.resize(nHalosBufferSend);// = (void *) malloc(locHalosBufferRecvSize);
//	locPartsBufferRecv = (void *) malloc(locPartsBufferRecvSize);
	barrMpi = MPI_Barrier(MPI_COMM_WORLD);

	MPI_Sendrecv(&locHalosBufferSend[0], locHalosBufferSend.size(), MPI_BYTE, recvTask, 0, 
		     &locHalosBufferRecv[0], locHalosBufferRecv.size(), MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &status);

//	MPI_Sendrecv(locPartsBufferSend, locPartsBufferSendSize, MPI_BYTE, recvTask, 0, 
//		     locPartsBufferRecv, locPartsBufferRecvSize, MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &status);

	cout << "On task=" << locTask << ", expecting =" << locHalosBufferRecvSize << " bytes from task=" << sendTask << endl;

	/* Now that the halos and particles have been communicated, clean the buffers */
	//free(locHalosBufferSend);
	//free(locPartsBufferSend);

	locHalosBufferNew.resize(nHalosBufferRecv);	
	locPartsBufferNew.resize(nHalosBufferRecv);	

	barrMpi = MPI_Barrier(MPI_COMM_WORLD);

	/* Unpack the buffers into the newly allocated Halos and Parts */
	for (int iH =0; iH < nHalosBufferRecv; iH++)
	{
		MPI_Unpack(&locHalosBufferRecv[0], locHalosBufferRecvSize, &bufferPosHaloRecv, &locHalosBufferNew[iH], 
	  		   sizeHalo, MPI_BYTE, MPI_COMM_WORLD);
		
		int nPartsThisHalo = locHalosBufferNew[iH].nPart;

		cout << "Task=" << locTask << " n=" << iH << " nPart=" <<  nPartsThisHalo << endl;

		//locPartsBufferNew[iH].resize(nPartsThisHalo);	

		//MPI_Unpack(locPartsBufferRecv, locPartsBufferRecvSize, &bufferPosPartRecv, &locPartsBufferNew[iH][0], 
		//	sizePart * nPartsThisHalo, MPI_BYTE, MPI_COMM_WORLD);

		locHalosBufferNew[iH].Info(); 
	}

	/* Now force the freeing of memory for all the recv buffers */
//	free(locHalosBufferRecv);
//	free(locPartsBufferRecv);
	
#endif		// Test Block
};



