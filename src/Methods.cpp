#include <algorithm>
#include <vector>

#include "Methods.h"
#include "general.h"

using namespace std;


Methods::Methods()
{
};


Methods::~Methods()
{
};

	
// Given two halos, decide whether to compare their particle content or not
bool Methods::CompareHalos(int iHalo, int jHalo)
{
	float rMax = 0.0;
	int iOne = 0, iTwo = 0;
	iOne = iUseCat; 
	iTwo = iUseCat + 1 % 1;

	rMax = locHalos[iOne][iHalo].rVir + locHalos[iTwo][jHalo].rVir; rMax *= dMaxFactor;

	// Only check for pairwise distance
	if (locHalos[iOne][iHalo].Distance(locHalos[iTwo][jHalo].X) < rMax)
		return true;
	else 
		return false;
	
};


void Methods::FindProgenitors()
{
	int nStepsCounter = floor(nLocHalos[iUseCat] / 50.);
	vector<int> nCommon;

	int iOne = 0, iTwo = 0;
	iOne = iUseCat; 
	iTwo = iUseCat + 1 % 1;

		for (int j = 0; j < nLocHalos[iOne]; j++)
		{
			if (j == nStepsCounter * floor(j / nStepsCounter))
				if (locTask == 0)			
					cout << "." << flush; 

			for (int i = 0; i < nLocHalos[iTwo]; i++)
				if (CompareHalos(j, i))
					nCommon = CommonParticles(locParts[iOne][j], locParts[iTwo][i]);
		}

};


vector<int> Methods::CommonParticles(vector<vector<unsigned long long int>> partsHaloOne, 
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

