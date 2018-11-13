
#include <math.h>

#include "Halo.h"
#include "utils.h"
#include "spline.h"
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
void Cosmology::GravAcc(Halo useHalo, float a0, float a1)
{
	float useM, deltaT, deltaV;
	float *rNorm;	rNorm = (float *) calloc(3, sizeof(float));

	/* Percentage of accuracy on the interpolated velocity computation */
	deltaV = VectorModule(useHalo.V) * 0.1;
	deltaT = A2Sec(a0, a1);

	/* It is faster to simply start the loop on all halos */
	for (int iH = 0; iH < locHalos[0].size(); iH++)
	{
		Halo thisHalo; thisHalo = locHalos[0][iH];
		float rHalos = useHalo.Distance(thisHalo.X);
		
		/* The cutoff radius depends on the mass of each halo */
		float rCutoff = sqrt((G_MpcMsunS * thisHalo.mTot * deltaT) / deltaV); 

		/* Now use the halo for comparison only if it is within a given radius */
		if (rHalos < rCutoff)
		{
			float thisAcc = G_MpcMsunS * thisHalo.mTot / (rHalos * rHalos);

			for (int iX = 0; iX < 3; iX++) 
				rNorm[iX] = (useHalo.X[iX] - thisHalo.X[iX]);	

			rNorm = UnitVector(rNorm);

			for (int iX = 0; iX < 3; iX++) 
				useHalo.V[iX] += rNorm[iX] * thisAcc;

			// TODO do some kind of leapfrog interpolation to compute the position at this intermediate step

		}
		
	}

};



float Cosmology::A2Sec(float a0, float a1)
{
	float time;

	// function a0, a1

	return time;
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
