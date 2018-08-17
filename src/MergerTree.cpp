#include <string>
#include <vector>
#include <algorithm>

#include "MergerTree.h"
#include "Halo.h"

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
	float rMax = 0.0;

	// do some check - if jHalo > nLocHalos ---> go look into the buffer halos FIXME

	rMax = locHalos[iOne][iHalo].rVir + locHalos[iTwo][jHalo].rVir; rMax *= dMaxFactor;

	// Only check for pairwise distance
	if (locHalos[iOne][iHalo].Distance(locHalos[iTwo][jHalo].X) < rMax)
		return true;
	else 
		return false;
	
};


/* Find all the progenitors (descendants) of the halos in catalog iOne (iTwo) */
void FindProgenitors(int iOne, int iTwo)
{
	int nStepsCounter = floor(nLocHalos[iUseCat] / 50.);
	vector<int> nCommon;

		for (int i = 0; i < nLocHalos[iOne]; i++)
		{
			if (i == nStepsCounter * floor(i / nStepsCounter) && locTask == 0)
					cout << "." << flush; 

			// FIXME instead of looping on all haloes ONLY select those within neighbouring nodes
			for (int j = 0; j < nLocHalos[iTwo]; j++)
				if (CompareHalos(i, j, iOne, iTwo))
					nCommon = CommonParticles(locParts[iOne][i], locParts[iTwo][j]);
		}
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

			//cout << "Task=" << locTask << " type=" << iT << " nPart=" << partsHaloOne[iT].size() << endl;
			//cout << "Task=" << locTask << " type=" << iT << " nPart=" << partsHaloTwo[iT].size() << endl;
	
			// Find out how many particles are shared among the two arrays
			iter = set_intersection(partsHaloOne[iT].begin(), partsHaloOne[iT].end(), 
				partsHaloTwo[iT].begin(), partsHaloTwo[iT].end(), thisCommon.begin());	

			// Resize the array and free some memory
			thisCommon.resize(iter - thisCommon.begin());
			//thisCommon.shrink_to_fit();		

			// Now compute how many particles in common are there
			nCommon[iT] = thisCommon.size();
		
			//cout << "Task=" << locTask << " type=" << iT << " nCommon=" << nCommon[iT] << endl;
	
			// Clear the vector and free all the allocated memory
			thisCommon.clear();
	 		thisCommon.shrink_to_fit();
		}
	}
	
	return nCommon;

};

