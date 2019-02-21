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
		GlobalGrid[iUseCat].AssignToGrid(locBuffHalos[iH].X, -iH-1);	// iH is negative - this is used for halos on the buffer, 
										// the -1 is added to avoid overlap with halo num 0

#ifdef VERBOSE
	cout << "Gathered " << locBuffHalos.size() << " halos in the buffer on task=" << locTask << endl;
	cout << "Gathered " << iBuffTotPart << " parts in the buffer on task=" << locTask << endl;
#else
	if (locTask == 0)
		cout << "Gathered halo and particle buffers. " << endl;
#endif
};


/* Once the merger trees are built, we gather the results on the master task.
 * This way we make sure that halos on the buffer scattered among multiple tasks 
 * are correctly syncronized and accounted for. */
void Communication::GatherMergerTrees(int iMTree)
{
	int nSendProgs = 0, nSendMains = 0, nRecvProgs= 0, nRecvMains = 0;
	int *sizeProgs = nullptr, *sizeMains = nullptr, *dispProgs = nullptr, *dispMains = nullptr;
	int *sizePTProgs = nullptr, *dispPTProgs = nullptr; 
	int *recvTrackProgs = nullptr, *recvTrackNComm = nullptr;
	Halo *recvMainHalos = nullptr, *recvProgHalos = nullptr;
	
	/* These are useful only non non-master task */
	vector<int> trackProgs, trackNComm;
	vector<Halo> mainHalos, progHalos;

	/* All tasks except the master now pack the halos for gathering */
	if (locTask != 0)
	{
		for (int iL = 0; iL < locMTrees[iMTree].size(); iL ++)
		{
			MergerTree thisTree = locMTrees[iMTree][iL];
			mainHalos.push_back(thisTree.mainHalo);
			trackProgs.push_back(thisTree.progHalo.size());			
	
			for (int iP = 0; iP < thisTree.progHalo.size(); iP++)
			{
				progHalos.push_back(thisTree.progHalo[iP]);
				
				for (int iT = 0; iT < nPTypes; iT++)
					trackNComm.push_back(thisTree.nCommon[iT][iP]);
			}
		}
		
		/* Count the objects to be sent and gathered */
		nSendProgs = progHalos.size();
		nSendMains = mainHalos.size();
	}

	/* First communicate how many main branch haloes are there and how many progenitors */
	MPI_Reduce(&nSendMains, &nRecvMains, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&nSendProgs, &nRecvProgs, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (locTask == 0)
	{
		cout << "Gathering " << nRecvProgs << " progenitors and " << nRecvMains << " main branch halos. " << endl;
		
		sizeProgs = (int *) malloc(totTask * sizeof(int));
		sizeMains = (int *) malloc(totTask * sizeof(int));
		dispProgs = (int *) malloc(totTask * sizeof(int));
		dispMains = (int *) malloc(totTask * sizeof(int));
		dispPTProgs = (int *) malloc(totTask * sizeof(int));
		sizePTProgs = (int *) malloc(totTask * sizeof(int));

		recvMainHalos = (Halo *) malloc(nRecvMains * sizeof(Halo));
		recvProgHalos = (Halo *) malloc(nRecvProgs * sizeof(Halo));
		recvTrackProgs = (int *) malloc(nRecvMains * sizeof(int));
		recvTrackNComm = (int *) malloc(nRecvProgs * nPTypes * sizeof(int));
	}

	/* Gather step 0 */
	MPI_Gather(&nSendMains, 1, MPI_INT, sizeMains, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(&nSendProgs, 1, MPI_INT, sizeProgs, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (locTask == 0)
	{
		nSendMains = 0;
		dispProgs[0] = 0;
		dispMains[0] = 0;
		sizePTProgs[0] = 0;
		dispPTProgs[0] = 0;

		for (int iT = 1; iT < totTask; iT++)
		{
			dispProgs[iT] = dispProgs[iT -1] + sizeProgs[iT -1];
			dispMains[iT] = dispMains[iT -1] + sizeMains[iT -1];

			/* Correct the sizes for number of particle types */
			dispPTProgs[iT] = dispProgs[iT] * nPTypes;
			sizePTProgs[iT] = sizeProgs[iT] * nPTypes;
		}
	}

	/* How many progenitors per main halo */
	MPI_Gatherv(&trackProgs[0], nSendMains, MPI_INT, recvTrackProgs, sizeMains, dispMains, MPI_INT, 0, MPI_COMM_WORLD);

	/* Particles shared by each progenitor */
	MPI_Gatherv(&trackProgs[0], nSendProgs * nPTypes, MPI_INT, recvTrackNComm, sizePTProgs, dispPTProgs, MPI_INT, 0, MPI_COMM_WORLD);

	if (locTask == 0)
	{
		for (int iT = 0; iT < totTask; iT++)
		{
			dispMains[iT] *= (int) sizeHalo;
			sizeMains[iT] *= (int) sizeHalo;
			dispProgs[iT] *= (int) sizeHalo;
			sizeProgs[iT] *= (int) sizeHalo;
		}
	
		mainHalos.clear();
		mainHalos.shrink_to_fit();
	}	

	/* Gather descendant main halos */
	MPI_Gatherv(&mainHalos[0], nSendMains * sizeHalo, MPI_BYTE, recvMainHalos, sizeMains, dispMains, MPI_BYTE, 0, MPI_COMM_WORLD);

	/* Gather progenitors */
	MPI_Gatherv(&progHalos[0], nSendProgs * sizeHalo, MPI_BYTE, recvProgHalos, sizeProgs, dispProgs, MPI_BYTE, 0, MPI_COMM_WORLD);

	/* Free the buffer and the bufffer information */
	if (locTask != 0)
	{
		mainHalos.clear();
		progHalos.clear();
	} else {
		free(sizeMains); free(dispMains);
		free(sizeProgs); free(dispProgs);
		free(sizePTProgs); free(dispPTProgs);
	}

	 /* Update the merger tree maps on the main task  */
	if (locTask == 0)
	{
		int iTrack = 0, iComm = 0, iAppend = 0;	

		for (int iMain = 0; iMain < nRecvMains; iMain++)
		{
			MergerTree thisTree;
			thisTree.mainHalo = recvMainHalos[iMain];

			int nProg = recvTrackProgs[iMain], iProg = 0;

			for (int iProg = 0; iProg < nProg; iProg++)
			{
				Halo thisProg = recvProgHalos[iTrack];
				thisTree.idProgenitor.push_back(thisProg.ID);
				thisTree.progHalo.push_back(thisProg);

				for (int iC = 0; iC < nPTypes; iC++)
				{
					thisTree.nCommon[iC].push_back(recvTrackNComm[iComm]);
					iComm++;
				}
			
				iTrack++;
			}
		
			/* if 0 --> There is no buffer exchange on 0 so just add everything to locMTree[0] */
			if (iMTree == 0)
			{
				/* After the 1 --> 0 FindProgenitor() comparison thisMap refers to 1 and nextMap to 0 */
				locMTrees[iMTree].push_back(thisTree);
				nextMapTrees[thisTree.mainHalo.ID] = locMTrees[iMTree].size()-1;

			/* Need to synchronize and update halos on the buffer */
			} else if (iMTree == 1) {
				map<uint64_t, int>::iterator iter;
				iter = thisMapTrees.find(thisTree.mainHalo.ID);
	
				/* This halo is already on the local buffer */
				if (iter != thisMapTrees.end())
				{
					int thisIndex = thisMapTrees[thisTree.mainHalo.ID];
					locMTrees[iMTree][thisIndex].Append(thisTree);		
					locMTrees[iMTree][thisIndex].SortByMerit();		
					iAppend++;
				} else {
					locMTrees[iMTree].push_back(thisTree);
					thisMapTrees[thisTree.mainHalo.ID] = locMTrees[1].size()-1;
				}
			}
		}	// for iMain 
	
		/* Free memory on task 0 */
		free(recvMainHalos); 
		free(recvProgHalos); 
		free(recvTrackProgs); 
		free(recvTrackNComm); 

	} 	// locTask == 0
}

/* Task 0 holds all the information about orphan halos. 
 * We MPI_Bcast it to all tasks so that each task can keep track of its orphans at the next step */
void Communication::SyncOrphanHalos()
{
	int nOrphIDs = 0, nLocOrphans = 0;

	if (locTask == 0)
		cout << "Synchronizing orphan halos across tasks ..." << endl;


	if (locTask == 0)
		nOrphIDs = allOrphIDs.size();	
	else
		nOrphIDs = 0;

	MPI_Bcast(&nOrphIDs, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (locTask != 0)
		allOrphIDs.resize(nOrphIDs);

	MPI_Bcast(&allOrphIDs[0], nOrphIDs * sizeof(uint64_t), MPI_BYTE, 0, MPI_COMM_WORLD);

	/* Each task looks for the orphan halo it holds */
	for (auto const& thisOrphID: allOrphIDs)
	{
		map<uint64_t, int>::iterator iter;

		iter = nextMapTrees.find(thisOrphID);
	
		/* If the ID is on this task, then copy the halo to the orphan list */
		if (iter != nextMapTrees.end())
		{
			int iTree = nextMapTrees[thisOrphID];
	
			/* Make sure we only track halos that were originally assigned to this task */
			if (iTree < nLocHalos[0])
			{
				Halo thisHalo = locHalos[0][iTree];
				thisHalo.isToken = true;
				thisHalo.nOrphanSteps++;

				locOrphHalos.push_back(thisHalo);

				/* Update the particle content */
				locOrphParts.push_back(locParts[0][iTree]);
				locOrphParts[nLocOrphans].resize(nPTypes);

				/* Copy particles by particle type */
				for (int iP = 0; iP < nPTypes; iP++)
					locOrphParts[nLocOrphans][iP] = locParts[0][iTree][iP];
					
			} // If nLocHalos[0]

			nLocOrphans++;	
		}
	}

	allOrphIDs.clear();
	allOrphIDs.shrink_to_fit();

	if (locTask == 0)
		cout << "Done." << endl;
}


/* 
 * Once the merger trees are built, we need to communicate the results 
 */
void Communication::SyncMergerTreeBuffer()
{
	/* Communicate the progenitor list for the halos in the buffer and keep track of their indexing */
	vector<Halo> buffSendProg, buffRecvProg, totBuffRecvProg;
	vector<int> buffSendProgIndex, buffRecvProgIndex, totBuffRecvProgIndex;
	vector<uint64_t> buffSendDescID, buffRecvDescID, totBuffRecvDescID;
	vector<int> buffSendComm, buffRecvComm, totBuffRecvComm;
	vector<int> trackHaloTask; 

	/* Tasks receiving and sending messages */
	int recvTask = 0, sendTask = 0;

	/* These buffers hold the number of halos to be sent/recvd to/from multiple tasks */
	int iBuffTotDescID = 0, nBuffSendDescID = 0, nBuffRecvDescID = 0, nTotBuffDescID = 0, nBuffRecvProg = 0, nBuffSendProg = 0;

	/* Determine the order of sending and receiving tasks to avoid gridlocks and make it consistent through
	 * all the tasks  */
	SetSendRecvTasks();

	if (locTask == 0)
		cout << "Synchronizing merger trees in the buffer regions..." << endl; 

	/* Loop on all the tasks except the local one */
	for (int iT = 0; iT < totTask-1; iT++)
	{
		sendTask = sendTasks[iT];
		recvTask = recvTasks[iT];

		/* At this step, locTask will send nBuffSendHalos to sendTask */
		int nBuffTotDesc = buffIndexSendHalo[sendTask].size();

#ifdef VERBOSE
		if (nBuffSendDescID == 0)
			cout << "Task=" << locTask << " has zero send buffer to " << sendTask << endl;  
#endif

		int iDesc = 0;

		/* Communicate the structure of these descendant halos (in backward mode) */
		for (int iP = 0; iP < nBuffTotDesc; iP ++)
		{
			int iH = buffIndexSendHalo[sendTask][iP];
			
			/* If the halo has no likely progenitor on this task then it's pointless to communicate it */
			if (locMTrees[1][iH].progHalo.size() > 0) 
			{
				buffSendDescID.push_back(locMTrees[1][iH].mainHalo.ID);
				buffSendProgIndex.push_back(locMTrees[1][iH].progHalo.size());

				for (int iProg = 0; iProg < buffSendProgIndex[iDesc]; iProg++)
				{
					Halo thisProg = locMTrees[1][iH].progHalo[iProg]; 
					buffSendProg.push_back(thisProg); 

					for (int iC = 0; iC < nPTypes; iC++)
						buffSendComm.push_back(locMTrees[1][iH].nCommon[iC][iProg]);
				}

				iDesc++;
			}
		}

		/* Size of all the progenitor halos to be sent */
		nBuffSendDescID = buffSendDescID.size();
		nBuffSendProg = buffSendProg.size();

#ifdef VERBOSE
		cout << iT << ") On task=" << locTask << ") sending " << nBuffSendProg << " halos to " << sendTask <<endl;
#endif

		/* Communicate the number of halo IDs on the buffer (progenitor only) */
		MPI_Sendrecv(&nBuffSendDescID, 1, MPI_INT, sendTask, 0, 
			     &nBuffRecvDescID, 1, MPI_INT, recvTask, 0, MPI_COMM_WORLD, &status);

		/* Communicate the number of progenitor halos */
		MPI_Sendrecv(&nBuffSendProg, 1, MPI_INT, sendTask, 0, 
			     &nBuffRecvProg, 1, MPI_INT, recvTask, 0, MPI_COMM_WORLD, &status);

		/* Multiply by two to take into account progenitor and descendant */
		size_t buffSendSizeDescID = nBuffSendDescID * sizeof(uint64_t);
		size_t buffRecvSizeDescID = nBuffRecvDescID * sizeof(uint64_t);

		buffRecvProgIndex.resize(nBuffRecvDescID);
		buffRecvDescID.resize(nBuffRecvDescID);	

#ifdef VERBOSE
		cout << iT << ") On task=" << locTask << ") recving " << nBuffRecvProg << " halos from " << recvTask <<endl;
		cout << iT << ") On task=" << locTask << ") recving " << buffRecvSizeDescID << " b from " << recvTask <<endl;
#endif
		MPI_Sendrecv(&buffSendDescID[0], buffSendSizeDescID, MPI_BYTE, sendTask, 0, 
			     &buffRecvDescID[0], buffRecvSizeDescID, MPI_BYTE, recvTask, 0, MPI_COMM_WORLD, &status);

		MPI_Sendrecv(&buffSendProgIndex[0], nBuffSendDescID, MPI_INT, sendTask, 0, 
			     &buffRecvProgIndex[0], nBuffRecvDescID, MPI_INT, recvTask, 0, MPI_COMM_WORLD, &status);

		buffRecvProg.resize(nBuffRecvProg);
		buffRecvComm.resize(nBuffRecvProg * nPTypes);

		/* Send the progenitor halos */
		MPI_Sendrecv(&buffSendProg[0], nBuffSendProg * sizeHalo, MPI_BYTE, sendTask, 0, 
			     &buffRecvProg[0], nBuffRecvProg * sizeHalo, MPI_BYTE, recvTask, 0, MPI_COMM_WORLD, &status);

		/* Send also the list of particle shared with each progenitor */
		MPI_Sendrecv(&buffSendComm[0], nBuffSendProg * nPTypes, MPI_INT, sendTask, 0, 
			     &buffRecvComm[0], nBuffRecvProg * nPTypes, MPI_INT, recvTask, 0, MPI_COMM_WORLD, &status);

		/* Append halo ids for the descendants and indexes for the progenitors to the total recv buffer */
		for (int iH = 0; iH < nBuffRecvDescID; iH++)
		{
			totBuffRecvDescID.push_back(buffRecvDescID[iH]);
			totBuffRecvProgIndex.push_back(buffRecvProgIndex[iH]);
		}

		int jH = 0;

		/* Append number of progenitors list */
		for (int iH = 0; iH < nBuffRecvProg; iH++)
		{
			totBuffRecvProg.push_back(buffRecvProg[iH]);

			for (int iC = 0; iC < nPTypes; iC++)
			{
				totBuffRecvComm.push_back(buffRecvComm[jH]);		
				jH++;
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		/* Clean the send buffer */
		if (buffRecvDescID.size() > 0)
		{
			buffSendProgIndex.clear();
			buffSendDescID.clear();
			buffSendProg.clear();
			buffSendComm.clear();
			buffSendProgIndex.shrink_to_fit();
			buffSendDescID.shrink_to_fit();
			buffSendProg.shrink_to_fit();
			buffSendComm.shrink_to_fit();
		}

		/* Clean the recv buffer */
		if (buffRecvProg.size() > 0)
		{
			buffRecvProgIndex.clear();
			buffRecvDescID.clear();
			buffRecvProg.clear();
			buffRecvComm.clear();
			buffRecvProgIndex.shrink_to_fit();
			buffRecvDescID.shrink_to_fit();
			buffRecvProg.shrink_to_fit();
			buffRecvComm.shrink_to_fit();
		}
	}	// Loop on the send/recv tasks

	int iProg = 0, iComm = 0;
	
	/* Now synchronize the connections of the halos on the buffer */
	for (int iH = 0; iH < totBuffRecvDescID.size(); iH++)
	{	
		uint64_t descID = totBuffRecvDescID[iH];
		int thisIndex = thisMapTrees[descID];
		int nProgs = totBuffRecvProgIndex[iH];

		for (int jH = 0; jH < nProgs; jH++)
		{

			Halo thisProg = totBuffRecvProg[iProg];
			locMTrees[1][thisIndex].idProgenitor.push_back(thisProg.ID);
			locMTrees[1][thisIndex].progHalo.push_back(thisProg);
	
			for (int iC = 0; iC < nPTypes; iC++)
			{
				locMTrees[1][thisIndex].nCommon[iC].push_back(totBuffRecvComm[iComm]);
				iComm++;
			}

			iProg++;
		}

		if (locMTrees[1][thisIndex].progHalo.size() > 1)
			locMTrees[1][thisIndex].SortByMerit();
	}

	/**/

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

