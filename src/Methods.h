/*
 * Wrapper class for different mathematical methods and algorithms used throughout the program
 */

#ifndef METHODS_H
#define METHODS_H

#include <vector>

#include "general.h" 

using namespace std;


class Methods {

public:
	Methods();
	~Methods();

	float dMaxFactor;
	bool fwdComparison;

	// Pairwise comparison of halos
	void FindProgenitors(int, int);

	// Given two (sorted) vectors, compare their content and return the number of common elements
	vector<int> CommonParticles(vector<vector<unsigned long long int>>, vector<vector<unsigned long long int>>);
	
	// Decide whether to compare two halos
	bool CompareHalos(int, int);

};
#endif
