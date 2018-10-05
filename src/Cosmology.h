#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include <iostream>
#include "Halo.h"
//#include "global_vars.h"


class Cosmology {

public:
	Cosmology();
	~Cosmology();

	void SetPlanck();
	void SetWMAP7();
	void SetArbitrary();	// TODO	

	
	float A2Sec(float, float);
	float Rho0(float, int);
	float RhoC(float, int);
		
	/* This function computes the gravitational acceleration due to the neighbouring halos for a token halo */
	void GravAcc(Halo, float, float);

	/* Cosmological parameters */
	float omegaL;
	float omegaM;
	float omegaDM;
	float omegaB;
	float h; 

	/* Useful units */
	float Mpc2km = 3.085e+19;
	float G_MpcMsunS = 4.517e-48; 
	float G_MKgS = 6.674e-11; 
	float Myr2s = 3.153e+13;
	float H100s = 3.241e-18;
	float Msun2g = 1.988e+33;

};
#endif
