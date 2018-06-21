#include <mpi.h>


class Communication {

public:
	// Empty constructor and destructor
	Communication() { };
	~Communication() { };	

	void ComputeBufferSize(void);
	void BufferSendRecv(void);


};

