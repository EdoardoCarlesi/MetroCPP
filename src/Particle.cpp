#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "Particle.h"

//using namespace std


Particle::Particle()
{
	X = new float[3];
}


Particle::~Particle()
{
	delete X;
}

