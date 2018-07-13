#ifndef COMMUNICATION_H
#define COMMUNICATION_H
#include <mpi.h>


class Communication {

public:
	// Empty constructor and destructor
	Communication() { };
	~Communication() { };	

	void BroadcastAndGatherGrid(void);
	void ComputeBufferSize(void);
	void BufferSendRecv(void);	

private:
	void ExchangeBuffers(void);
	
	void SetSendRecvTasks(void);

	// Buffers for halo and node communication
	vector<vector<int>> buffIndexNodeHalo;
	
	// Here we store the full list of halo indexes
	vector<vector<int>> buffIndexSendHalo;

	// This only holds the number of haloes to be received
	vector<int> sendTasks;
	vector<int> recvTasks;
};
#endif
