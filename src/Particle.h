#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>
#include "general.h"

class Particle {

public:
	Particle();
	~Particle();

	short int Type;	
	unsigned long long int ID;

#ifdef PART_SNAPSHOT
	float *X;
	float *V;
#endif
};


#endif
