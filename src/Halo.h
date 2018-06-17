#ifndef HALO_H
#define HALO_H

#include "Particle.h"

class Halo {

public:
	Halo();
	~Halo();

	float Distance(float *);
	float *getPos(void) { return X; };
	float *getVel(void) { return V; };

	
private:
	float *X;
	float *V;
	float M;
	float R;

	int NPart;
	long int ID;

	Particle *Part;
};

#endif
