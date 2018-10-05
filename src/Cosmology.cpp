
#include <math.h>

#include "utils.h"
#include "global_vars.h"
#include "Cosmology.h"


Cosmology::Cosmology()
{
	/* Default constructor sets Planck parameters */
	SetPlanck();	
};


Cosmology::~Cosmology()
{};



float Cosmology::Rho0(float boxSizeMpc, int nPart)
{
	float fact0 = 100. / 256.;
	float mass0 = 1.05217e+11 / 20.;
	float fact1 = boxSizeMpc / nPart;
	float mass1 = pow(fact1/fact0, 3) * mass0;
	float rho0 = mass1 * nPart / pow(boxSizeMpc, 3);

	return rho0;
};



/* Compute the gravitational acceleration around a given halo */
float Cosmology::GravAcc()
{
};



float Cosmology::RhoC(float boxSizeMpc, int nPart)
{
	return Rho0(boxSizeMpc, nPart) * (1.0 / (1.0 - omegaL));
};



void Cosmology::SetPlanck()
{
	omegaDM = 0.26;
	omegaL = 0.69;
	omegaM = 0.31;
	omegaB = 0.05;
	h = 0.67;
};


	
void Cosmology::SetWMAP7()
{
	omegaDM = 0.23;
	omegaL = 0.73;
	omegaM = 0.27;
	omegaB = 0.04;
	h = 0.7;
};
