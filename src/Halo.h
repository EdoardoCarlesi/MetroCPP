#ifndef HALO_H
#define HALO_H

#include <string>
#include <cstddef>
#include <iostream>
#include <malloc.h>
#include "Cosmology.h"

using namespace std;


class Halo {

public:
	/* Standard constructor / destructor */
	Halo();
	~Halo();

	float mTot, mDM, mGas, mFM, mStar, rVir;
	float cNFW, rsNFW;			// NFW parameters
	float lambda, vMax, sigV;
	float X[3], V[3], L[3];
	float fMhires;

	// Token halos keeping track of "lost" subhalos need to set this to TRUE
	bool isToken;
	
	// Number of time steps during which the halo has been "orphan" of a progenitor
	int nOrphan;
	
	// Number of subhalos
	int nSub;
	
	// Total number of particles is set to 7
	int nPart[NPTYPES+1];
	unsigned long long int ID, hostID;

	// Compute the halo distance from a given point
	float Distance(float *);
	
	/* NFW Mass at a given radius */
	float M_NFW(float);
	
	// Velocity (module) of the halo wrt a given velocity vector
	float RelativeVelocity(float *);

	// Cout some infos on the halo (velocity, position, etc.)
	void Info(void);
};

#endif
