#include <mpi.h>


class Communication {

public:
	// Empty constructor and destructor
	Communication() { };
	~Communication() { };	

	void BroadcastAndGatherGrid(void);
	void ComputeBufferSize(void);
	void BufferSendRecv(void);
};

