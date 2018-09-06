#include <string>
#include <vector>
#include <algorithm>

#include "MergerTree.h"
#include "Halo.h"

#include "utils.h"
#include "global_vars.h"

using namespace std;



MergerTree::MergerTree()
{
};


MergerTree::~MergerTree()
{
};


void MergerTree::sortByMerit(void)
{
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
	vector<int> nCommon, indexes;
	int totCmp = 0, totSkip = 0, totNCommon = 0;

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

			for (int j = 0; j < locHalos[iTwo].size(); j++)
			{
				nCommon = CommonParticles(locParts[iOne][thisIndex], locParts[iTwo][j]);
				//cout << j << " onTask= " << locTask << ", npart= " << locParts[iOne][thisIndex][1].size() << endl;
			}
	
			if (locTask == 1)
				cout << i << " " << thisIndex << " " << nCommon[1] << endl;
					
			if (nCommon[1] > 10) 
			{		
				totNCommon += nCommon[1];
				totCmp++;
			}

	}

#else	/* Standard comparison */

		for (int i = 0; i < nLocHalos[iOne]; i++)
		{
			if (i == nStepsCounter * floor(i / nStepsCounter) && locTask == 0)
					cout << "." << flush; 
#ifdef CMP_ALL	/* Compare ALL the halos located on the task - used only as a benchmark */

			for (int j = 0; j < locHalos[iTwo].size(); j++)
				nCommon = CommonParticles(locParts[iOne][i], locParts[iTwo][j]);
						
			if (nCommon[1] > 10) 
			{		
				totNCommon += nCommon[1];
				totCmp++;
			}

			for (int j = 0; j < locBuffHalos.size(); j++)
				nCommon = CommonParticles(locParts[iOne][i], locBuffParts[-j]);
						
			if (nCommon[1] > 10) 
			{		
				totNCommon += nCommon[1];
				totCmp++;
			}

#else		/* Compare to a subset of halos */

			Halo thisHalo = locHalos[iOne][i];
			rSearch = facRSearch * thisHalo.rVir;

			/* We only loop on a subset of halos */
			indexes = GlobalGrid[iTwo].ListNearbyHalos(thisHalo.X, rSearch);

			for (int j = 0; j < indexes.size(); j++)
			{
				int k = indexes[j];

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
						totCmp++;
						totNCommon += nCommon[1];
					} else {
						totSkip++;
					}
				} // Halo Comparison

			}	// for j, k = index(j)

#endif		// compare all halos
		} // for i halo, the main one

#endif		// ifdef ZOOM


		if (locTask == 0)
			cout << "\n" << locTask << ") TotComp: " << totCmp << ", totComm: " << totNCommon << " skip: " << totSkip << endl; 
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

