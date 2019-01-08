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

#ifndef COMMUNICATION_H
#define COMMUNICATION_H
#include <mpi.h>


class Communication {

public:
	// Empty constructor and destructor
	Communication() { };
	~Communication() { };	

#ifdef ZOOM
	void SyncOrphanHalos(void);
#else 
	void SyncIndex(void);
	void BroadcastAndGatherGrid(void);
	void SyncMergerTreeBuffer(void);
#endif

	void BufferSendRecv(void);	
	
	void CleanBuffer(void);

private:
	// Determine the send and recv tasks consistently
	vector<int> sendTasks;
	vector<int> recvTasks;

#ifndef ZOOM
	void SetSendRecvTasks(void);

	// Communicate the buffers across all tasks
	void ExchangeBuffers(void);
#endif
	// Buffers for halo and node communication
	vector<vector<int>> buffIndexNodeHalo;
	
	// Here we store the full list of halo indexes
	vector<vector<int>> buffIndexSendHalo;
};
#endif
