#include <math.h>
#include <iostream>

using namespace std;

#include "Particle.h"
#include "Halo.h"


Halo::Halo()
{
	Mtot = 0.0; 	Mgas = 0.0; 	Mdm = 0.0; 	NPart = 0;
	Rvir = 0.0;	Spin = 0.0;
	ID = 0;
};


Halo::~Halo()
{
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


void Halo::Info(void)
{
	cout << "Task " << locTask << " Halo ID " << ID << endl;
	printf("Mdm: %.3e, Mgas: %.3e, Npart: %d\n", Mdm, Mgas, NPart);

};

