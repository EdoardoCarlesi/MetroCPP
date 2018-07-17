#ifndef MERGERTREE_H
#define MERGERTREE_H

#include <string>

#include "Halo.h"

using namespace std;



//class Halo;


class MergerTree {

public:

	MergerTree();
	~MergerTree();

	int *nProgenitors; // At each step store the number of progenitors found

	Halo *mainHalo;		// At each step stores one main halo
	Halo **haloProgenitors;	// Multiple array containing at each step all progenitors

	void sortByMerit(void);	// Once possible progenitors have been found, compare

};


/*
 * Functions used to build the merger trees 
 */

// Pairwise comparison of halos
void FindProgenitors(int, int);

// Given two (sorted) vectors, compare their content and return the number of common elements
vector<int> CommonParticles(vector<vector<unsigned long long int>>, vector<vector<unsigned long long int>>);
	
// Decide whether to compare two halos
bool CompareHalos(int, int, int, int);





#endif
