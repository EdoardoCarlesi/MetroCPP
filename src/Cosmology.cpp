
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
void Cosmology::GravAcc(int haloIndex, float a0, float a1)
{
	float rHalo[3], rHalos, rCutoff, *rNorm;	
	float a_i0[3], a_i1[3], x_i0[3];
	float useM, deltaT, deltaV;

#ifdef TEST
	/* Percentage of accuracy on the interpolated velocity computation */
	deltaV = VectorModule(useHalo.V) * 0.1;
	deltaT = A2Sec(a0, a1);

	/* It is faster to simply start the loop on all halos */
	for (int iH = 0; iH < locCleanTrees[iNumCat-1].size(); iH++)
	{
		Halo thisHaloA1; thisHaloA1 = locCleanTrees[iNumCat-1][iH];
		rHalos = useHalo.Distance(thisHaloA1.X);
		
		/* The cutoff radius depends on the mass of each halo */
		rCutoff = sqrt((G_MpcMsunS * thisHaloA1.mTot * deltaT) / deltaV); 

		/* Now use the halo for comparison only if it is within a given radius */
		if (rHalos < rCutoff)
		{
			thisAcc = G_MpcMsunS * thisHaloA1.mTot / (rHalos * rHalos);

			for (int iX = 0; iX < 3; iX++) 
				rHalo[iX] = (useHalo.X[iX] - thisHalo.X[iX]);	

			rNorm = UnitVector(&rHalo[0]);

			for (int iX = 0; iX < 3; iX++) 
				a_i1[iX] += rNorm[iX] * thisAcc;

		}

		/* For the leapfrog algorithm we need to run the code at the two timesteps */
		Halo thisHaloA1; thisHaloA1 = locHalos[1][iH];
		rHalos = useHalo.Distance(thisHaloA1.X);
		
		/* The cutoff radius depends on the mass of each halo */
		rCutoff = sqrt((G_MpcMsunS * thisHaloA1.mTot * deltaT) / deltaV); 

		/* Now use the halo for comparison only if it is within a given radius */
		if (rHalos < rCutoff)
		{
			//float thisAcc = G_MpcMsunS * thisHalo.mTot / (rHalos * rHalos);

			for (int iX = 0; iX < 3; iX++) 
				rHalo[iX] = (useHalo.X[iX] - thisHalo.X[iX]);	

			rNorm = UnitVector(&rHalo[0]);

			for (int iX = 0; iX < 3; iX++) 
				a_i0[iX] += rNorm[iX] * thisAcc;

		}

	
	}
#endif
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
