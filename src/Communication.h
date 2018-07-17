#ifndef COMMUNICATION_H
#define COMMUNICATION_H
#include <mpi.h>


class Communication {

public:
	// Empty constructor and destructor
	Communication() { };
	~Communication() { };	

	void BroadcastAndGatherGrid(void);
	void BufferSendRecv(void);	

private:
	// Determine the send and recv tasks consistently
	void SetSendRecvTasks(void);
	vector<int> sendTasks;
	vector<int> recvTasks;

	// Communicate the buffers across all tasks
	void ExchangeBuffers(void);

	// Buffers for halo and node communication
	vector<vector<int>> buffIndexNodeHalo;
	
	// Here we store the full list of halo indexes
	vector<vector<int>> buffIndexSendHalo;
};
#endif
