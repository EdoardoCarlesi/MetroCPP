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
	for (int iM = 0; iM < mTree.size(); iM++)
		mTree[iM].Clean();

	mTree.clear();
	mainHalo.clear();
};


void MergerTree::Info()
{
	//if (nPart > 1000)
	{
		//cout << "Task=" << locTask << " " << idDescendant << " " << idProgenitor.size() << endl;
		cout << "Task=" << locTask << " " << idDescendant << " " << idProgenitor.size() << endl;

		if (tokenProgenitor)
			for (int iP = 0; iP < idProgenitor.size(); iP++)
				cout << "Ind=" << indexProgenitor[iP] << " ID= " << idProgenitor[iP] << endl;
		else
			for (int iP = 0; iP < idProgenitor.size(); iP++)
				cout << "Ind=" << indexProgenitor[iP] << " ID= " << idProgenitor[iP] << " NP=" << nCommon[1][iP] << endl;

		//cout << "--" << endl;
	}
};


void MergerTree::Clean()
{
	nCommon.clear();
	idProgenitor.clear();
	indexProgenitor.clear();
};


MergerTree::MergerTree()
{
	//tokenProgenitor = false;
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

	cmpHalo = locHalos[iTwo][jHalo];
	rMax = locHalos[iOne][iHalo].rVir + cmpHalo.rVir;
	rMax *= 20.0; 
	//rMax = locHalos[iOne][iHalo].rVir + cmpHalo.rVir; 
	//rMax *= dMaxFactor * 5.0;
	
	//cout << iHalo << " " << jHalo << " " << locHalos[iOne][iHalo].Distance(cmpHalo.X) << " " << rMax << endl;

	// Only check for pairwise distance
	//if (locHalos[iOne][iHalo].Distance(cmpHalo.X) < rMax)
	//	return true;
	//else 
	//	return false;
#else
	// TODO introduce more checks on the time, velocity and so on
	// do some check - if jHalo > nLocHalos ---> go look into the buffer halos FIXME
	if (jHalo >= 0)
		cmpHalo = locHalos[iTwo][jHalo];
	if (jHalo < 0)
		cmpHalo = locBuffHalos[-jHalo];

	rMax = locHalos[iOne][iHalo].rVir + cmpHalo.rVir; 
#endif

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

	locMTrees[iOne].clear();
	locMTrees[iOne].shrink_to_fit();
	locMTrees[iOne].resize(nLocHalos[iOne]);

	//cout << locTask << ") nHalos= " << nLocHalos[iOne] << endl;

#ifdef ZOOM
	vector<int> locHaloPos;
	int thisIndex = 0;
	int halosPerTask = 0, halosRemaind = 0;

	if (iOne < iTwo)	// Only forward comparison is split, when doing backward comparison each task need to 
				// control all the halos, otherwise the trees will not be built correctly (in ZOOM mode)
	{
		halosPerTask = int(nLocHalos[iOne] / totTask);
		halosRemaind = nLocHalos[iOne] % totTask;

		/* The first halosRemaind tasks get one halo more */
		if (locTask < halosRemaind)
			halosPerTask += 1;

		cout << "Task=" << locTask << " is looping on " << halosPerTask << " halos in the forward loop." << endl;
	} else {		// If doing backward correlation we have to loop on all halos (ZOOM mode)
		halosPerTask = locTreeIndex.size();
		
		/*
		cout << "Task=" << locTask << " is looping on " << halosPerTask << " halos in the backward loop." << endl;
		cout <<  locHalos[iOne].size() << ", nHalos=" << nLocHalos[iOne] << ", indexMax=" << locTreeIndex[halosPerTask-1] << endl;
		cout <<  locHalos[iTwo].size() << ", nHalos=" << nLocHalos[iTwo] << ", indexMax=" << locTreeIndex[halosPerTask-1] << endl;
		*/
	}

	for (int i = 0; i < halosPerTask; i++)
	{
		if (iOne < iTwo)
			thisIndex = locTask + i * totTask;
		else
			thisIndex = locTreeIndex[i];
			//thisIndex = i;

		Halo thisHalo = locHalos[iOne][thisIndex];
		
		/* Save some halo properties on the mtree vector */
		locMTrees[iOne][thisIndex].idDescendant = thisHalo.ID;
		locMTrees[iOne][thisIndex].nPart = thisHalo.nPart[nPTypes];
		locMTrees[iOne][thisIndex].nCommon.resize(nPTypes);

			/* In a zoom-in run, we loop over all the halos on the iTwo step */
			for (int j = 0; j < locHalos[iTwo].size(); j++)
			{ 
                                //if (CompareHalos(i, j, iOne, iTwo))
				{
					nCommon = CommonParticles(locParts[iOne][thisIndex], locParts[iTwo][j]);
	
					int totComm = 0;

					for (int iC = 0; iC < nCommon.size(); iC++)
						totComm += nCommon[iC];

					/* This is very important: we keep track of the merging history ONLY if the number 
					   of common particles is above a given threshold */
					if (totComm > minPartCmp) 
					{		
						for (int iT = 0; iT < nPTypes; iT++)
						{
							totNCommon[iT] += nCommon[iT];
							locMTrees[iOne][thisIndex].nCommon[iT].push_back(nCommon[iT]);
						}

						locMTrees[iOne][thisIndex].idProgenitor.push_back(locHalos[iTwo][j].ID);
						locMTrees[iOne][thisIndex].indexProgenitor.push_back(j);
					
						if (locMTrees[iOne][thisIndex].idDescendant == locHalos[iTwo][j].ID)
						{
							locMTrees[iOne][thisIndex].tokenProgenitor == true;
							//cout << "Found orphan progenitor " << locHalos[iTwo][j].ID << endl;
						}

						/* We keep track of the halos on iTwo that have been matched on iOne on the 
						 * local task, so that we can avoid looping on all iTwo halos afterwards */
						if (iOne < iTwo)
							locTreeIndex.push_back(j);
					
						totCmp++;
					}
				} // CompareHalos
			}

			/*
			if (locHalos[iOne][thisIndex].isToken && iOne < iTwo && locMTrees[iOne][thisIndex].idProgenitor.size() > 0)
			{
				cout << "**** on Task= " << locTask << " found a progenitor for orphan halo: " << endl; 
				locHalos[iOne][thisIndex].Info();
				locHalos[iTwo][locMTrees[iOne][thisIndex].indexProgenitor[0]].Info();
				locMTrees[iOne][thisIndex].Info();
				cout << "****" << endl;
			}
			*/

			/* Very important: if it turns out the halo has no likely progenitor, and has a number of particles above 
			 * minPartHalo, then we keep track of its position in the global locHalos[iOne] array
			 */
			if (locMTrees[iOne][thisIndex].idProgenitor.size() == 0 && 
				locMTrees[iOne][thisIndex].nPart > minPartHalo && iOne < iTwo)
				{
					for (int iT = 0; iT < nPTypes; iT++)
						locMTrees[iOne][thisIndex].nCommon[iT].push_back(locHalos[iOne][thisIndex].nPart[iT]);

				//locTreeIndex.push_back(thisIndex);
					orphanHaloIndex.push_back(thisIndex);
					locMTrees[iOne][thisIndex].tokenProgenitor = true;
					//locMTrees[iOne][thisIndex].idProgenitor.push_back(locHalos[iOne][thisIndex].ID);
					//locMTrees[iOne][thisIndex].indexProgenitor.push_back();
				} else {
					locMTrees[iOne][thisIndex].tokenProgenitor = false;
				}
			/*
			if (locMTrees[iOne][thisIndex].idProgenitor.size() == 0 && 
				locMTrees[iOne][thisIndex].nPart > minPartHalo && iOne < iTwo)
			{
				cout << "\n----Orphan halo on task=" << locTask << ", index= " << thisIndex << endl;
				locHalos[iOne][thisIndex].Info();
				cout << "----" << endl;
			}	
			*/

#ifdef DEBUG		// Sanity check
			if (locMTrees[iOne][thisIndex].nPart > 1000)
				locMTrees[iOne][thisIndex].Info();
#endif
	}	// End loop on iOne halo

	/* Sort and clean of multiple entries the locTreeIndex */
	sort(locTreeIndex.begin(), locTreeIndex.end());
	locTreeIndex.erase(unique(locTreeIndex.begin(), locTreeIndex.end()), locTreeIndex.end());

#else	/* Standard comparison */

		for (int i = 0; i < nLocHalos[iOne]; i++)
		{
			locMTrees[iOne][i].idDescendant = locHalos[iOne].ID;
			locMTrees[iOne][i].nPart = locHalos[iOne].nPart[nPTypes];
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

				locMTrees[iOne][i].idDescendant = locHalos[iOne].ID;
				locMTrees[iOne][i].nPart = locHalos[iOne].nPart[nPTypes];
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

/*
		if (locTask == 0)
			cout << "\n" "Compared a total of : " << totCmp 
				<< " halos, with a total number of: " << totNCommon[1] << " Type 1 DM particles in common and "
				<< "  a total number of: " << totNCommon[2] << " Type 2 DM particles in common. " << endl; 
*/
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



void InitTrees(int nUseCat)
{
	locCleanTrees.resize(nUseCat-1);
	allHalos.resize(nUseCat);
	copy(locHalos[0].begin(), locHalos[0].end(), back_inserter(allHalos[iNumCat]));
};



void CleanTrees(int iStep)
{
	// TODO implement the sort by merit function!!!!!
	for (int iSim = 0; iSim < 2; iSim++)
		for (int iTree = 0; iTree < locMTrees[iSim].size(); iTree++)
			locMTrees[iSim][iTree].sortByMerit();

	if (locTask == 0)
		cout << "Cleaning Merger Tree connections, back and forth..." << endl;

	for (int iTree = 0; iTree < locMTrees[0].size(); iTree++)
	{
		unsigned long long int mainID = locHalos[0][iTree].ID;

		MergerTree mergerTree;
		mergerTree.nCommon.resize(nPTypes);
		mergerTree.idDescendant = mainID;
		mergerTree.tokenProgenitor = locMTrees[0][iTree].tokenProgenitor;
		mergerTree.nPart = locHalos[0][iTree].nPart[nPTypes];
	
		/* At each step we only record the connections between halos in catalog 0 and catalog 1, without attempting at a
		   reconstruction of the full merger history. This will be done later. */
		for (int iProg = 0; iProg < locMTrees[0][iTree].idProgenitor.size(); iProg++)
		{
			int jTree = locMTrees[0][iTree].indexProgenitor[iProg];
			unsigned long long int progID = locMTrees[0][iTree].idProgenitor[iProg];
			unsigned long long int descID = locMTrees[1][jTree].idProgenitor[0];	// The progenitor in the backwards tree

#ifdef DEBUG
			// Useful for debugging
			if (descID == 0)
			{
				cout << "WARNING. Progenitor ID not assigned: " << progID << " " << descID 
					<< " | " << iTree << " " << jTree << endl;

				locHalos[0][iTree].Info();
				locHalos[1][jTree].Info();
			}
#endif

			if (mainID == descID)
			{
				if (locMTrees[0][iTree].tokenProgenitor)
				{
					mergerTree.idProgenitor.push_back(locMTrees[0][iTree].idDescendant);
					//mergerTree.indexProgenitor.push_back(locMTrees[0][iTree].indexProgenitor[0]);

					for(int iT = 0; iT < nPTypes; iT++)
						mergerTree.nCommon[iT].push_back(locMTrees[0][iTree].nCommon[iT][iProg]);

				} else {
					mergerTree.idProgenitor.push_back(progID);
					mergerTree.indexProgenitor.push_back(jTree);

					for(int iT = 0; iT < nPTypes; iT++)
						mergerTree.nCommon[iT].push_back(locMTrees[0][iTree].nCommon[iT][iProg]);
				}
			}
		}

/*
		if (mergerTree.tokenProgenitor)
		{
			//locMTrees[0][iTree].Info();
			//cout << locMTrees[0][iTree].idDescendant << endl;
			//cout << locMTrees[0][iTree].idProgenitor.size() << endl;
			//mergerTree.idProgenitor[0] = locMTrees[0][iTree].idDescendant;
			//cout << "token " << mergerTree.idProgenitor.size() << endl;
		}
*/

		if (mergerTree.idProgenitor.size() > 0)
			locCleanTrees[iStep-1].push_back(mergerTree);

		mergerTree.Clean();
	}
	
	locTreeIndex.clear();
	locTreeIndex.shrink_to_fit();
};



void DebugTrees()
{
	if (locTask == 0)
		cout << "Debugging trees for steps=" << locCleanTrees.size() << endl;
 
	for (int iC = 0; iC < locCleanTrees.size(); iC++)
	{	
		cout << "Task=(" << locTask << ") Size=(" << locCleanTrees[iC].size() << ") Step=(" << iC << ") " << endl;

		for (int iT = 0; iT < locCleanTrees[iC].size(); iT++)
				locCleanTrees[iC][iT].Info();
	}
};



