#include <math.h>

#include "Halo.h"


float Halo::Distance(float *Pos)
{
	int iX = 0;
	float dX = 0.0;

	for (iX = 0; iX < 3; iX++)
	{
		dX += pow(Pos[iX] - Halo::X[iX], 2);	
	}

return dX;

};

