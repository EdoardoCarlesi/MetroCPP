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


#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include <iostream>
#include "Halo.h"
#include "spline.h"
#include "global_vars.h"

//using namespace tk;


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
	
	float H2t(float);
	
	/* This function computes the gravitational acceleration due to the neighbouring halos for a token halo */
	void GravAcc(int, float, float);

	/* Cosmological parameters */
	float omegaL;
	float omegaM;
	float omegaDM;
	float omegaB;
	float h; 

	/* Splines used for interpolation */
	tk::spline pk;
	tk::spline a;

private:

	float Kick();
	float Drift();

	float InitH2t();

	/* Useful units */
	double Mpc2km = 3.085e+19;
	double G_MpcMsunS = 4.517e-48; 
	double G_MKgS = 6.674e-11; 
	double Myr2s = 3.153e+13;
	double H100s = 3.241e-18;
	double Msun2g = 1.988e+33;

};
#endif
