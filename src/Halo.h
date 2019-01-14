/*
 *   METROC++: MErger TRees On C++, a scalable code for the computation of merger trees in cosmological simulations.
 *   Copyright (C) Edoardo Carlesi 2018-2019
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


#ifndef HALO_H
#define HALO_H

#include <string>
#include <cstddef>
#include <iostream>
#include <malloc.h>
//#include "Cosmology.h"

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
	int nOrphanSteps;
	
	// Number of subhalos
	int nSub;
	
	// Total number of particles is set to normal types +1
	int nPart[NPTYPES+1];

	uint64_t ID, hostID;

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
