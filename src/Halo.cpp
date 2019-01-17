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


/*
 * Halo.cpp defines the basic Halo class, with some very simple properties.
 * For the ease of communication, it is better to avoid dynamically allocated data, so that it is easier to 
 * MPI_Sendrecv Halo-type arrays in the buffers.
 */

#include <math.h>
#include <malloc.h>
#include <iostream>
#include <string>

#include "Cosmology.h"
#include "global_vars.h"
#include "Halo.h"

using namespace std;


Halo::Halo()
{
	mTot = 0.0; 	mGas = 0.0; 	mDM = 0.0; 	
	rVir = 0.0;	lambda = 0.0;   nOrphanSteps = 0;
	ID = 0;		hostID = 0;
	fMhires = 1.000;
	isToken = false; 	
	
	nPart[nPTypes] = 0;

	for (int iX = 0; iX <3; iX++)
	{
		X[iX] = -1.0;
		V[iX] = 0.00;
	}

	for (int iT = 0; iT < nPTypes; iT++)
		nPart[0] = 0;
};


Halo::~Halo()
{
	//delete Part;
};


float Halo::Distance(float *Pos)
{
	int iX = 0;
	float dX = 0.0;

	for (iX = 0; iX < 3; iX++)
	{
		dX += pow(Pos[iX] - X[iX], 2);	
	}
	
	dX = sqrt(dX);
	
	return dX;
};


float Halo::M_NFW(float r)
{	
	//float = 4.0 * 3.14159 // Cosmology --> ;
};


void Halo::Info(void)
{
	if (isToken == 0)
		cout << "Task " << locTask << " (TOKEN) Halo ID " << ID << endl;
	else
		cout << "Task " << locTask << " Halo ID " << ID << endl;
		
	printf("Mtot: %.3e, ID: %lu, Npart: %d, X:(%.2f, %.2f, %.2f), V:(%.2f, %.2f, %.2f), fMhires:%.3f, isToken:%d\n", 
		mTot, ID, nPart[nPTypes], X[0], X[1], X[2], V[0], V[1], V[2], fMhires, isToken);
};
