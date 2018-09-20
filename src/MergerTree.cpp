#include <string>
#include <vector>
#include <algorithm>
#include <math.h>

#include "MergerTree.h"
#include "Halo.h"

#include "utils.h"
#include "global_vars.h"

using namespace std;



HaloTree::HaloTree()
{
};


HaloTree::~HaloTree()
{
	Clean();
};


void HaloTree::Clean()
{
};



MergerTree::MergerTree()
{
};


MergerTree::~MergerTree()
{
};


void MergerTree::sortByMerit()
{
	// merit = pow( common/(nP1 * nP2), 2.0 )
};


       /************************************************************************* 
	* These are general functions and are not part of the class Merger Tree *
	*************************************************************************/
		
/* Given two halos, decide whether to compare their particle content or not */
bool CompareHalos(int iHalo, int jHalo, int iOne, int iTwo)
{
	float min = 100000.0, max = 0.0;
	float rMax = 0.0, vMax = 0.0, fVel = 0.5e-2, vOne = 0.0, vTwo = 0.0;
	Halo cmpHalo;

#ifdef ZOOM

#else
	// TODO introduce more checks on the time, velocity and so on
	// do some check - if jHalo > nLocHalos ---> go look into the buffer halos FIXME
	if (jHalo >= 0)
		cmpHalo = locHalos[iTwo][jHalo];
	if (jHalo < 0)
		cmpHalo = locBuffHalos[-jHalo];
#endif

	rMax = locHalos[iOne][iHalo].rVir + cmpHalo.rVir; 
	vOne = VectorModule(locHalos[iOne][iHalo].V); vTwo = VectorModule(cmpHalo.V);
	vMax = (vOne + vTwo) * fVel;
	rMax *= dMaxFactor * vMax;

	// Only check for pairwise distance
	if (locHalos[iOne][iHalo].Distance(cmpHalo.X) < rMax)
		return true;
	else 
		return false;
	
};


/* Find all the progenitors (descendants) of the halos in catalog iOne (iTwo) */
void FindProgenitors(int iOne, int iTwo)
{
	int nStepsCounter = floor(nLocHalos[iUseCat] / 50.);
	float rSearch = 0, facRSearch = 20.0;
	vector<int> nCommon, indexes, totNCommon;
	int totCmp = 0, totSkip = 0; 

	totNCommon.resize(nPTypes);
	locMTrees[iOne].resize(nLocHalos[iOne]);

	//cout << locTask << ") nHalos= " << nLocHalos[iOne] << endl;

#ifdef ZOOM
	vector<int> locHaloPos;
	int thisIndex = 0;
	int halosPerTask = 0, halosRemaind = 0;

	halosPerTask = int(nLocHalos[iOne] / totTask);
	halosRemaind = nLocHalos[iOne] % totTask;

	/* The first halosRemaind tasks get one halo more */
	if (locTask < halosRemaind)
		halosPerTask += 1;

	for (int i = 0; i < halosPerTask; i++)
	{
		thisIndex = locTask + i * totTask;
		Halo thisHalo = locHalos[iOne][thisIndex];

		locMTrees[iOne][thisIndex].nCommon.resize(nPTypes);

			for (int j = 0; j < locHalos[iTwo].size(); j++)
			{
				nCommon = CommonParticles(locParts[iOne][thisIndex], locParts[iTwo][j]);
	
				/* This is very important: we keep track of the merging history ONLY based on the number 
				   of common DM particles */
				if (nCommon[1] > 10) 
				{		
					for (int iT = 0; iT < nPTypes; iT++)
					{
						totNCommon[iT] += nCommon[iT];
						locMTrees[iOne][thisIndex].nCommon[iT].push_back(nCommon[iT]);
					}

					locMTrees[iOne][thisIndex].idProgenitor.push_back(locHalos[iTwo][j].ID);
					locMTrees[iOne][thisIndex].indexProgenitor.push_back(j);

					totCmp++;
				}

			}
	}

#else	/* Standard comparison */

		for (int i = 0; i < nLocHalos[iOne]; i++)
		{
			locMTrees[iOne][i].nCommon.resize(nPTypes);

			if (i == nStepsCounter * floor(i / nStepsCounter) && locTask == 0)
					cout << "." << flush; 
#ifdef CMP_ALL	/* Compare ALL the halos located on the task - used only as a benchmark */

			for (int j = 0; j < locHalos[iTwo].size(); j++)
			{
				nCommon = CommonParticles(locParts[iOne][i], locParts[iTwo][j]);
						
				if (nCommon[1] > 10) 
				{		
					for (int iT = 0; iT < nPTypes; iT++)
					{
						totNCommon[iT] += nCommon[iT];
						locMTrees[iOne][i].nCommon[iT].push_back(nCommon[iT]);
					}

					locMTrees[iOne][i].idProgenitor.push_back(locHalos[iTwo][j].ID);
					locMTrees[iOne][i].indexProgenitor.push_back(j);
					totCmp++;
				}
			}

			for (int j = 0; j < locBuffHalos.size(); j++)
			{
				nCommon = CommonParticles(locParts[iOne][i], locBuffParts[j]);
							
				if (nCommon[1] > 10) 
				{		
					for (int iT = 0; iT < nPTypes; iT++)
					{
						totNCommon[iT] += nCommon[iT];
						locMTrees[iOne][i].nCommon[iT].push_back(nCommon[iT]);
					}

					locMTrees[iOne][i].idProgenitor.push_back(locBuffHalos[j].ID);
					locMTrees[iOne][i].indexProgenitor.push_back(-j);
				}

				totCmp++;
			}

#else		/* Compare to a subset of halos */

			Halo thisHalo = locHalos[iOne][i];
			rSearch = facRSearch * thisHalo.rVir;

			/* We only loop on a subset of halos */
			indexes = GlobalGrid[iTwo].ListNearbyHalos(thisHalo.X, rSearch);

			/* The vector "indexes" contains the list of haloes (in the local memory & buffer) to be compared */
			for (int j = 0; j < indexes.size(); j++)
			{
				int k = indexes[j];

				locMTrees[iOne][i].nCommon.resize(nPTypes);

				/* Compare halos --> this functions checks whether the two halos are too far 
				   or velocities are oriented on opposite directions */
				if (CompareHalos(i, k, iOne, iTwo))
				{	
					if (k >= 0)
						nCommon = CommonParticles(locParts[iOne][i], locParts[iTwo][k]);
					else
						nCommon = CommonParticles(locParts[iOne][i], locBuffParts[-k]);

					if (nCommon[1] > 10) 
					{		
						for (int iT = 0; iT < nPTypes; iT++)
						{
							locMTrees[iOne][i].nCommon[iT].push_back(nCommon[iT]);
							totNCommon[iT] += nCommon[iT];
						}
	
					if (k >= 0)
					
						locMTrees[iOne][i].idProgenitor.push_back(locHalos[k].ID);
					else 
						locMTrees[iOne][i].idProgenitor.push_back(locBuffHalos[-k].ID);
						
					locMTrees[iOne][i].indexProgenitor.push_back(k);

						totCmp++;
					} else {
						totSkip++;
					}
				} // Halo Comparison

			}	// for j, k = index(j)

#endif		// compare all halos
		} // for i halo, the main one

#endif		// ifdef ZOOM


		if (locTask == 0)
			cout << "\n" << locTask << ") TotComp: " << totCmp << ", DM TotComm: " << totNCommon[1] << endl; 
};


/* Given a pair of haloes, determine the number of common particles */
vector<int> CommonParticles(vector<vector<unsigned long long int>> partsHaloOne, 
	vector<vector<unsigned long long int>> partsHaloTwo)
{
	vector<int> nCommon; 
	vector<unsigned long long int>::iterator iter;
	vector<unsigned long long int> thisCommon;

	nCommon.resize(nPTypes);	

	for (int iT = 0; iT < nPTypes; iT++)
	{
		int thisSize = partsHaloOne[iT].size();

		if (thisSize > 0)
		{
			//cout << locTask << " " << iT << " " << partsHaloOne[iT][0] << " " << partsHaloOne[iT][1]
			//	<< " " << partsHaloTwo[iT][0] << " " << partsHaloTwo[iT][1] << endl;

			// This is the maximum possible number of common particles
			thisCommon.resize(thisSize);
	
			// Find out how many particles are shared among the two arrays
			iter = set_intersection(partsHaloOne[iT].begin(), partsHaloOne[iT].end(), 
				partsHaloTwo[iT].begin(), partsHaloTwo[iT].end(), thisCommon.begin());	

			// Resize the array and free some memory
			thisCommon.resize(iter - thisCommon.begin());
			//thisCommon.shrink_to_fit();		

			// Now compute how many particles in common are there
			nCommon[iT] = thisCommon.size();
		
			// Clear the vector and free all the allocated memory
			thisCommon.clear();
	 		thisCommon.shrink_to_fit();
		}
	}
	
	return nCommon;

};



void CleanTrees(int iStep)
{
	if (iStep == 0)
		allHTrees.resize(1);

	for (int iSim = 0; iSim < 2; iSim++)
		for (int iTree = 0; iTree < locMTrees[iSim].size(); iTree++)
			locMTrees[iSim][iTree].sortByMerit();

	for (int iTree = 0; iTree < locMTrees[0].size(); iTree++)
	{
		unsigned long long int mainID = locHalos[0][iTree].ID;

		/* At each step we only record the connections between halos in catalog 0 and catalog 1, without attempting at a
		   reconstruction of the full merger history. This will be done later and saved into the allHaloTrees	*/
		HaloTree haloTree;
		haloTree.mTree.resize(1);
		haloTree.mainHalo.resize(1);
	
		haloTree.mTree[0].nCommon.resize(nPTypes);
		haloTree.mainHalo[0] = locHalos[0][iTree];

		for (int iProg = 0; iProg < locMTrees[0][iTree].indexProgenitor.size(); iProg++)
		{
			int jTree = locMTrees[0][iTree].indexProgenitor[iProg];
			unsigned long long int progID = locMTrees[1][jTree].idProgenitor[0];
		
			if (mainID == progID)
			{
				haloTree.mTree[0].idProgenitor.push_back(progID);	
				
				for(int iT = 0; iT < nPTypes; iT++)
					haloTree.mTree[0].nCommon[iT].push_back(locMTrees[0][iTree].nCommon[iT][iProg]);
			}

			locHTrees[iStep].push_back(haloTree);
		}

		haloTree.Clean();
	}

};


