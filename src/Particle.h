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
	float X[3];
	float V[3];
#endif
};


#endif
