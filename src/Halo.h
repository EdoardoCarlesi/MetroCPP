#ifndef HALO_H
#define HALO_H

#include "Particle.h"

class Halo {

public:
	Halo();
	~Halo();
	
	float X[3];
	float V[3];
	float M;
	float R;

	int NPart;
	long int ID;

	Particle *P;

	float Distance(float *);

};

#endif
