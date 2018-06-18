#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "Particle.h"

using namespace std;

Particle::Particle()
{
#ifdef PART_SNAPSHOT
	X = new float[3];
#endif
}


Particle::~Particle()
{
#ifdef PART_SNAPSHOT
	delete X;
#endif
}

