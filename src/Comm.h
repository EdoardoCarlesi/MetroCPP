#include <mpi.h>


class Comm {

public:
	// Empty constructor and destructor
	Comm() { };
	~Comm() { };	

	void ComputeBufferSize(void);
	void BufferSendRecv(void);


};

