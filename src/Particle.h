#ifndef PARTICLES_H
#define PARTICLES_H
#include <cmath>

class Particle {

public:
	Particle();
	~Particle();

	float *getPos(void) { return X; };
	int getType(void) { return Type; };

private:
	float *X;
	//float V[3];
	short int Type;	
	long int ID;
};


#endif
