#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "Particle.h"

using namespace std;



Particle::Particle()
{
	ID = 0; 	Type = 9; 
};


Particle::~Particle()
{
};


void Particle::ReadLineAHF(const char * lineRead)
{
	        sscanf(lineRead, "%llu %hu", &ID, &Type); 
};
