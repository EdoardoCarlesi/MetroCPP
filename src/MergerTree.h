#ifndef MERGERTREE_H
#define MERGERTREE_H

#include <string>

#include "Halo.h"
#include "Particle.h"

using namespace std;



class Halo;


class MergerTree {

public:

	MergerTree();
	~MergerTree();

	int *nProgenitors; // At each step store the number of progenitors found

	Halo *mainHalo;		// At each step stores one main halo
	Halo **haloProgenitors;	// Multiple array containing at each step all progenitors

	void sortByMerit(void);	// Once possible progenitors have been found, compare

};





#endif
