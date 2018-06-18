#ifndef HALO_H
#define HALO_H

#include <string>
#include <cstddef>
#include <iostream>

//#include "general.h"
#include "Particle.h"

using namespace std;

class Particle;

class Halo {

public:
	// Standard constructor / destructor
	Halo();
	~Halo();

	float X[3];
	float V[3];
	float L[3];

	float Mtot;
	float Mgas;
	float Mdm;
	float Mstar;
	float Rvir;

	float Spin;
	float Vmax;

	unsigned int NPart;
	unsigned long long int ID;

	Particle *Part;

	// Compute the halo distance from a given point
	float Distance(float *);
		
	// Velocity (module) of the halo wrt a given velocity vector
	float RelativeVelocity(float *);

	// Set methods: sets the halo properties reading from different formats of halo catalogs
	void readLineAHF(string);
	void readLineFOF(string);

	void Info(void);
};

#endif
