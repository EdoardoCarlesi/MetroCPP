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

	locHalosBufferSendSize = nHalosBufferSend * sizeHalo;
	locHalosBufferSend.resize(nHalosBufferSend);

	locPartsBufferSend.resize(nHalosBufferSend);
	locPartsBufferSendSize = 0;

	cout << "Particle task=" << locTask << ", h1 " << locParts[0][0].ID << endl;	
	cout << "Particle task=" << locTask << ", h2 " << locParts[1][0].ID << endl;	
	cout << "Particle task=" << locTask << ", h3 " << locParts[2][0].ID << endl;	

	for (int iP = 0; iP < nHalosBufferSend; iP ++)
	{
		locHalosBufferSend[iP] = locHalos[iP];
		locPartsBufferSendSize += locHalos[iP].nPart * sizePart;
		locPartsBufferSend[iP].resize(locHalos[iP].nPart);
		locPartsBufferSend[iP] = locParts[iP];
	}

	// FIXME this is just a temporary setting
	/* RecvTask is recieving from LocTask and SendTask is sending to LocTask */
	recvTask = (locTask + totTask - 1) % (totTask);
	sendTask = (locTask + 1) % (totTask);

	//cout << "On task=" << locTask << " sending " << locHalosBufferSendSize << " bytes of halos to task=" << recvTask << endl;
	//cout << "On task=" << locTask << " sending " << locHalosBufferSend.size() * sizeHalo << " bytes to task=" << recvTask << endl;
 
	cout << "On task=" << locTask << " sending " << locPartsBufferSendSize << " bytes of parts to task=" << recvTask << endl;
	cout << "On task=" << locTask << " sending " << locPartsBufferSend.size() * sizePart << " bytes to task=" << recvTask << endl;
 
	cout << "On task=" << locTask << " ";
	locHalos[0].Info();

	/* Communicate halo and particle numbers to be sent across tasks */
	MPI_Sendrecv(&nHalosBufferSend, 1, MPI_INT, recvTask, 0, &nHalosBufferRecv, 1, MPI_INT, sendTask, 0, MPI_COMM_WORLD, &status);

	MPI_Sendrecv(&locHalosBufferSendSize, sizeof(size_t), MPI_BYTE, recvTask, 0, 
		     &locHalosBufferRecvSize, sizeof(size_t), MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &status);

	//cout << "On task=" << locTask << ", expecting=" << nHalosBufferRecv << " halos from task=" << sendTask << endl;
	//cout << "On task=" << locTask << ", expecting=" << locHalosBufferRecvSize << " bytes from task=" << sendTask << endl;

	locHalosBufferRecv.resize(nHalosBufferRecv);

	barrMpi = MPI_Barrier(MPI_COMM_WORLD);

	MPI_Sendrecv(&locHalosBufferSend[0], locHalosBufferSendSize, MPI_BYTE, recvTask, 0, 
		     &locHalosBufferRecv[0], locHalosBufferRecvSize, MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &status);

	//cout << "AFTER SENDRECV On task=" << locTask << " ";
	//locHalosBufferRecv[0].Info();

	locPartsBufferRecv.resize(nHalosBufferRecv);

	for (int iP = 0; iP < nHalosBufferRecv; iP ++)
		locPartsBufferRecv[iP].resize(locHalosBufferRecv[iP].nPart);

	MPI_Sendrecv(&locPartsBufferSendSize, sizeof(size_t), MPI_BYTE, recvTask, 0, 
		     &locPartsBufferRecvSize, sizeof(size_t), MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &status);

	MPI_Sendrecv(&locPartsBufferSend[0][0], locPartsBufferSendSize, MPI_BYTE, recvTask, 0, 
		     &locPartsBufferRecv[0][0], locPartsBufferRecvSize, MPI_BYTE, sendTask, 0, MPI_COMM_WORLD, &status);

	cout << "Particle task=" << locTask << ", sending =" << locPartsBufferSendSize << " MB to task=" << recvTask << endl;
	cout << "Particle task=" << locTask << ", expecting=" << locPartsBufferRecvSize << " MB from task=" << sendTask << endl;

	cout << "H-Buffer on task=" << locTask << ", " << locHalosBufferSend.size() << endl;
	cout << "P-Buffer on task=" << locTask << ", " << locPartsBufferSend[0].size() << endl;

#ifdef TEST_BLOCK	// FIXME clean up the memory...

	/* Now that the halos and particles have been communicated, clean the buffers */
	locHalosBufferSend.clear();
	locHalosBufferSend.shrink_to_fit();

	//cout << "H-Buffer on task=" << locTask << ", " << locHalosBufferSend.size() << endl;
	//cout << "P-Buffer on task=" << locTask << ", " << locPartsBufferSend[0].size() << endl;

	for (int iP = 0; iP < nHalosBufferSend; iP ++)	
	{
		locPartsBufferSend[iP].clear();
		locPartsBufferSend[iP].shrink_to_fit();
	}

	locPartsBufferSend.clear();
	locPartsBufferSend.shrink_to_fit();
	
#endif		// Test Block
};



