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
	Communication() {};
	~Communication();	

#ifndef ZOOM
	void BroadcastAndGatherGrid(void);
	void SyncMergerTreeBuffer(void);
	void GatherMergerTrees(int);
#endif

	void BufferSendRecv(void);	
	
	void CleanBuffer(void);

private:
	/* Store the send and recv tasks consistently */
	vector<int> sendTasks;
	vector<int> recvTasks;

#ifndef ZOOM
	void SetSendRecvTasks(void);

	void ExchangeBuffers(void);
#endif
	/* Which nodes lie on which tasks, and which halo they contain */
	vector<vector<int>> buffIndexNodeHalo;
	vector<vector<int>> buffIndexSendHalo;
};
#endif
