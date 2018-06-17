#ifndef PARTICLES_H
#define PARTICLES_H
#include <cmath>

class Particle {

public:
	Particle();
	~Particle();

	float X[3];
	//float V[3];
	short int Type;	
	long int ID;
};


#endif
