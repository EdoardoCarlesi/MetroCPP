#ifndef HALO_H
#define HALO_H

#include <string>
#include <cstddef>
#include <iostream>
#include <malloc.h>

#include "general.h"

using namespace std;

#define NPTYPES 6	/* This needs to be defined at compile time so that the Halo class has a well-defined size,
			 * which makes MPI communication easier
			 */


class Halo {

public:
	// Standard constructor / destructor
	Halo();
	~Halo();

	float mTot, mDM, mGas, mFM, mStar, rVir;
	float lambda, vMax, sigV;
	float X[3], V[3], L[3];

	// Token halos keeping track of "lost" subhalos need to set this to TRUE
	bool isToken;

	int nSub;
	
	// Total number of particles is set to 7
	int nPart[NPTYPES+1];
	unsigned long long int ID, hostID;

	// Compute the halo distance from a given point
	float Distance(float *);
		
	// Velocity (module) of the halo wrt a given velocity vector
	float RelativeVelocity(float *);

	// Set methods: sets the halo properties reading from different formats of halo catalogs
	void ReadLineAHF(const char *);
	
	// TODO Not implemented yet
	void ReadLineFOF(const char *);

	// Cout some infos on the halo (velocity, position, etc.)
	void Info(void);
};

#endif
