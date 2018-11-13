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
	if (isToken)
		cout << "Task " << locTask << " (TOKEN) Halo ID " << ID << endl;
	else
		cout << "Task " << locTask << " Halo ID " << ID << endl;
		
	printf("Mtot: %.3e, ID: %.llu, Npart: %d, X:(%.2f, %.2f, %.2f), V:(%.2f, %.2f, %.2f), fMhires:%.3f\n", 
		mTot, ID, nPart[nPTypes], X[0], X[1], X[2], V[0], V[1], V[2], fMhires);

};
