#ifndef MERGERTREE_H
#define MERGERTREE_H

#include <string>

#include "Halo.h"

using namespace std;

/* Due to the reliance of this class on vector template, we cannot directly MPI_Sendrecv the MergerTrees */
class MergerTree {

public:
	MergerTree();
	~MergerTree();

	int nPart;					// Number of particles per halo

	unsigned long long int idDescendant;		// ID of the descendant halo
	bool tokenProgenitor;				// If no progenitor is found, then create a token halo
							// with the same particle content to keep tracking it at subsequent steps
	vector<vector<int>> nCommon;			// Particles in common are separated per particle type
	vector<unsigned long long int> idProgenitor;	// IDs of progenitors
	vector<unsigned long long int> indexProgenitor;	// local array index of progenitors

	void sortByMerit(void);				// Once possible progenitors have been found, compare
	void Clean(void);
	void Info(void);
};


/* This class stores the main halo and its progenitors at each step */
class HaloTree {
	
public:
	HaloTree();
	~HaloTree();

	int nStep;
	
	vector<MergerTree> mTree; 
	vector<Halo> mainHalo;				// This traces the main progenitor

	void SmoothTree(void);				// Smooths over fly-bys 
	void FixTree(void);				// Looks for missing subhalos and fixes with token halos at the missing positions
	void Clean(void);

	void WriteMergerTree(void);			// Prints all the informations 
	void WriteTrajectory(void);
	void WriteMAH(void);
	void WriteIDs(void);				// Only print the ID of each halo and its main progenitor
};


/*
 * Functions used to build the merger trees 
 */

void CleanTrees(int);
void DebugTrees(void);

// Pairwise comparison of halos
void FindProgenitors(int, int);

// Decide whether to compare two halos
bool CompareHalos(int, int, int, int);

// Given two (sorted) vectors, compare their content and return the number of common elements
vector<int> CommonParticles(vector<vector<unsigned long long int>>, vector<vector<unsigned long long int>>);

#endif
