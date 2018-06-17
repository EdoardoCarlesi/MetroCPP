#include <math.h>

#include "Halo.h"


Halo::Halo()
{
	X = new float[3];
	V = new float[3];
};


Halo::~Halo()
{
	delete X;
	delete V;
	delete Part;
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

