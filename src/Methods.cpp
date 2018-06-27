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


vector<int> Methods::CommonParticles(vector<vector<unsigned long long int>> partsHaloOne, 
	vector<vector<unsigned long long int>> partsHaloTwo)
{
	vector<int> nCommon; 
	vector<unsigned long long int>::iterator iter;
	vector<unsigned long long int> thisCommon;
	int nPTypes = 0;

	nPTypes = partsHaloOne.size();
	nCommon.resize(nPTypes);	

	for (int iT = 0; iT < nPTypes; iT++)
	{
		int thisSize = partsHaloOne[iT].size();
		
		if (thisSize > 0)
		{
			thisCommon.resize(thisSize);

			cout << "Task=" << locTask << " type=" << iT << " nPart=" << partsHaloOne[iT].size() << endl;
			cout << "Task=" << locTask << " type=" << iT << " nPart=" << partsHaloTwo[iT].size() << endl;

			iter = set_intersection(partsHaloOne[iT].begin(), partsHaloOne[iT].end(), 
				partsHaloTwo[iT].begin(), partsHaloTwo[iT].end(), thisCommon.begin());	

			thisCommon.resize(iter - thisCommon.begin());
		
			// How many particles in common there are
			nCommon[iT] = thisCommon.size();
		
			cout << "Task=" << locTask << " type=" << iT << " nCommon=" << nCommon[iT] << endl;
	
			// Clear the vector
			thisCommon.clear();
	 		thisCommon.shrink_to_fit();
		}
	}
#ifdef TEST_BLOCK
#endif

};

