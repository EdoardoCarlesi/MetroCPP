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
	cout << "Task=" << locTask << " " << mainHalo.ID << " " << idProgenitor.size() << endl;

	if (isOrphan)
		for (int iP = 0; iP < idProgenitor.size(); iP++)
			cout << "Ind=" << indexProgenitor[iP] << " ID= " << idProgenitor[iP] << endl;
		else
			for (int iP = 0; iP < idProgenitor.size(); iP++)
				cout << "Ind=" << indexProgenitor[iP] << " ID= " << idProgenitor[iP] << " NP=" << nCommon[1][iP] << endl;

};



void MergerTree::Clean()
{
	for (int iC = 0; iC < nCommon.size(); iC++)
	{
		nCommon[iC].clear();
	}

	nCommon.clear();

	idProgenitor.clear();
	indexProgenitor.clear();
	subHalos.clear();
};



MergerTree::MergerTree()
{
};



MergerTree::~MergerTree()
{
};


/* This sorting algorithm might be very inefficient but it's straightforward to implement, 
 * plus we will rarely deal with halos with more than 10^6 progenitors */
void MergerTree::sortByMerit(int jSimu)
{
	vector<unsigned long long int> tmpIdx;
	vector<float> allMerit;
	vector<int> idx, tmpIndex;
	//vector<Halo> tmpSubHalos;

	int iSimu = (jSimu + 1) % 2; 
	float merit = 0.0;

	for (int iM = 0; iM < idProgenitor.size(); iM++)
	{
		int thisIndex = indexProgenitor[iM];
		Halo thisProgHalo = locHalos[iSimu][thisIndex];
		int nComm = 0;

		for(int iC = 0; iC < nPTypes; iC++)
				nComm += nCommon[iC][iM];

		merit = ((float) nComm * nComm) / (mainHalo.nPart[1] * thisProgHalo.nPart[1]);
		merit *= merit;

		//cout << iM << " " << merit << " " << thisProgHalo.mTot << " " << endl;

		allMerit.push_back(merit);
	}
	
	idx = SortIndexes(allMerit);
	tmpIndex.resize(idx.size());
	tmpIdx.resize(idx.size());
	//tmpSubHalos.resize(idx.size());

	for (int iM = 0; iM < idProgenitor.size(); iM++)
	{
		//cout << iM << ", " << idx[iM] << ", " << idProgenitor[iM] << ", " << idProgenitor.size() << endl;
		tmpIdx[iM] = idProgenitor[idx[iM]];
		tmpIndex[iM] = indexProgenitor[idx[iM]];
		//tmpSubHalos[iM] = subHalos[idx[iM]];
	}	

	idProgenitor = tmpIdx;
	indexProgenitor = tmpIndex;
	//subHalos = tmpSubHalos;

	/*	
	cout << idx << ", " << locHalos[iSimu][indexProgenitor[idx]].mTot << endl;
	cout << 0   << ", " << locHalos[iSimu][indexProgenitor[0]].mTot << endl;
	cout << "........." << endl;

	int idx = 0;
	for (int iM = 0; iM < idProgenitor.size(); iM++)
	{
		cout << locTask << "] (" << idx[iM] << ", " << iM << ") (" 
			<< allMerit[idx[iM]] << ", " << allMerit[iM] << ") (" 
			<< locHalos[iSimu][indexProgenitor[iM]].mTot << ", " 
			<< locHalos[iSimu][indexProgenitor[idx[iM]]].mTot << ") " << endl; 
	}
		cout << endl;
	}*/



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

#ifdef ZOOM	// TODO this is not clear yet
	cmpHalo = locHalos[iTwo][jHalo];
	rMax = locHalos[iOne][iHalo].rVir + cmpHalo.rVir;
	rMax *= 25.0; 
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

	// If we are dealing with token halos, enlarge the search radius
	if (locHalos[iOne][iHalo].isToken || cmpHalo.isToken)
			rMax *= 2.0;
	
	// Only check for pairwise distance
	if (locHalos[iOne][iHalo].Distance(cmpHalo.X) < rMax)
	{
		//cout << "Halo comparison: " << locHalos[iOne][iHalo].Distance(cmpHalo.X) << " r: " << rMax << endl;
		//locHalos[iOne][iHalo].Info();
		//cmpHalo.Info();
		//cout << "--" << endl;
		return true;
	} else {
		//cout << "Skipping halo comparison: " << endl;
		//locHalos[iOne][iHalo].Info();
		//cmpHalo.Info();
		//cout << "--" << endl;
		return false;
	}
	
};


/* Find all the progenitors (descendants) of the halos in catalog iOne (iTwo) */
void FindProgenitors(int iOne, int iTwo)
{
	vector<int> thisNCommon, indexes, totNCommon;
	float rSearch = 0, facRSearch = 20.0;
	int nStepsCounter = floor(nLocHalos[iUseCat] / 50.);
	int totCmp = 0, totSkip = 0, nLocOrphans = 0; 

	locMTrees[iOne].clear();
	locMTrees[iOne].shrink_to_fit();
	locMTrees[iOne].resize(nLocHalos[iOne]);

#ifdef ZOOM
	vector<int> locHaloPos;
	int thisIndex = 0, halosPerTask = 0, halosRemaind = 0;

	/* For consistency across tasks and to keep track of all possible orphans, IDs need to be initialized here */
	for (int iT = 0; iT < locHalos[iOne].size(); iT++)
		locMTrees[iOne][iT].mainHalo.ID = locHalos[iOne][iT].ID;

	/* The way forward/backward comparisons are done is different: in the forward mode, halos are split among the tasks as 
	 * Task 1 --> halo 1, task 2 ---> halo 2, assuming that the halos are ordered with increasing (decreasing) number of particles.
	 * In the backward mode, we assign to each task only those halos which have been matched in the fwd step, to avoid unnecessary
	 * comparisons.   */
	if (iOne < iTwo)	
	{
		halosPerTask = int(nLocHalos[iOne] / totTask);
		halosRemaind = nLocHalos[iOne] % totTask;
		locTreeIndex.clear();
		locTreeIndex.shrink_to_fit();

		/* The first halosRemaind tasks get one halo more */
		if (locTask < halosRemaind)
			halosPerTask += 1;

		if (locTask == 0)
			cout << "Looping on " << halosPerTask << " halos in the forward loop." << endl;
	} else {

		/* If doing backward correlation we have to loop on all halos (ZOOM mode) */
		halosPerTask = locTreeIndex.size();
		
		if (locTask == 0)
			cout << "Looping on " << halosPerTask << " halos in the backward loop." << endl;
	}

	for (int i = 0; i < halosPerTask; i++)
	{
		/* This distribution assumes that halos are ordered per number of particles, to improve the load balancing */
		if (iOne < iTwo)
			thisIndex = locTask + i * totTask;
		else
			thisIndex = locTreeIndex[i];

		Halo thisHalo = locHalos[iOne][thisIndex];
		
		/* Save some halo properties on the mtree vector */
		locMTrees[iOne][thisIndex].mainHalo = thisHalo;
		locMTrees[iOne][thisIndex].nCommon.resize(nPTypes);
	
			/* In a zoom-in run, we loop over all the halos on the iTwo step */
			for (int j = 0; j < locHalos[iTwo].size(); j++)
			{ 
				thisNCommon.resize(nPTypes);
                                //if (CompareHalos(i, j, iOne, iTwo))	// TODO this does not really improve speed, need to do more tests
				{
					thisNCommon = CommonParticles(locParts[iOne][thisIndex], locParts[iTwo][j]);
	
					int totComm = 0;
					totComm = thisNCommon[0] + thisNCommon[1] + thisNCommon[2];

					/* This is very important: we keep track of the merging history ONLY if the number 
					 * of common particles is above a given threshold */
					if (totComm > minPartCmp) 
					{	
						/* Common particles are separated per particle type */	
						for (int iT = 0; iT < nPTypes; iT++)
						{	
							//cout << j << ", " << iT << ", " << thisNCommon[iT] << endl;
							locMTrees[iOne][thisIndex].nCommon[iT].push_back(thisNCommon[iT]);
						}

						locMTrees[iOne][thisIndex].idProgenitor.push_back(locHalos[iTwo][j].ID);
						locMTrees[iOne][thisIndex].indexProgenitor.push_back(j);
					
						/* If the two halos have the same ID, we are dealing with a token halo
						 * placed to trace the progenitor of an orphan halo */
						if (locMTrees[iOne][thisIndex].mainHalo.ID == locHalos[iTwo][j].ID)
							locMTrees[iOne][thisIndex].isOrphan == true;

						/* We keep track of the halos on iTwo that have been matched on iOne on the 
						 * local task, so that we can avoid looping on all iTwo halos afterwards */
						if (iOne < iTwo)
							locTreeIndex.push_back(j);
					
						totCmp++;
				
					}

				} // CompareHalos
				
				thisNCommon.clear();
				thisNCommon.shrink_to_fit();
			}

			locMTrees[iOne][thisIndex].sortByMerit(iOne);

			/* Very important: if it turns out the halo has no likely progenitor, and has a number of particles above 
			 * minPartHalo, then we keep track of its position in the global locHalos[iOne] array  */
			if (locMTrees[iOne][thisIndex].idProgenitor.size() == 0 && 
				locMTrees[iOne][thisIndex].mainHalo.nPart[1] > minPartHalo && iOne < iTwo)
				{
					for (int iT = 0; iT < nPTypes; iT++)
						locMTrees[iOne][thisIndex].nCommon[iT].push_back(locHalos[iOne][thisIndex].nPart[iT]);

					orphanHaloIndex.push_back(thisIndex);
					locMTrees[iOne][thisIndex].isOrphan = true;
				} else {
					locMTrees[iOne][thisIndex].isOrphan = false;
				}

#ifdef DEBUG		// Sanity check
			if (locMTrees[iOne][thisIndex].nPart > 1000)
				locMTrees[iOne][thisIndex].Info();
#endif
	}	// End loop on iOne halo

	/* Sort and clean of multiple entries the locTreeIndex */
	sort(locTreeIndex.begin(), locTreeIndex.end());
	locTreeIndex.erase(unique(locTreeIndex.begin(), locTreeIndex.end()), locTreeIndex.end());



			/**************************************************
			 *	Comparison for fullbox simulations        *
			 **************************************************/
#else						



		for (int i = 0; i < nLocHalos[iOne]; i++)
		{
			locMTrees[iOne][i].mainHalo = locHalos[iOne][i];
			locMTrees[iOne][i].nCommon.resize(nPTypes);

			if (i == nStepsCounter * floor(i / nStepsCounter) && locTask == 0)
					cout << "." << flush; 

			Halo thisHalo = locHalos[iOne][i];
			rSearch = facRSearch * thisHalo.rVir;

			/* We only loop on a subset of halos */
			indexes = GlobalGrid[iTwo].ListNearbyHalos(thisHalo.X, rSearch);

		// TODO: this does not take into account (yet) the orphan halos in FULLBOX simulations
		// Now, in FullBox mode, Orphan Halos do not need to be communicated through all the tasks,
		// Which makes things easier to synchronize
	
			/* The vector "indexes" contains the list of haloes (in the local memory & buffer) to be compared */
			for (int j = 0; j < indexes.size(); j++)
			{
				int k = indexes[j];

				//locMTrees[iOne][i].idDescendant = locHalos[iOne].ID;
				//locMTrees[iOne][i].nPart = locHalos[iOne].nPart[nPTypes];
				//locMTrees[iOne][i].mainHalo = locHalos[iOne];
				//locMTrees[iOne][i].nCommon.resize(nPTypes);

				/* Compare halos --> this functions checks whether the two halos are too far 
				   or velocities are oriented on opposite directions */
				if (CompareHalos(i, k, iOne, iTwo))
				{	
					if (k >= 0)
						thisNCommon = CommonParticles(locParts[iOne][i], locParts[iTwo][k]);
					else
						thisNCommon = CommonParticles(locParts[iOne][i], locBuffParts[-k]);

					int totComm = 0;
					totComm = thisNCommon[0] + thisNCommon[1] + thisNCommon[2];

					/* This is very important: we keep track of the merging history ONLY if the number 
					 * of common particles is above a given threshold */
					if (totComm > minPartCmp) 
					{		
						for (int iT = 0; iT < nPTypes; iT++)
						{
							locMTrees[iOne][i].nCommon[iT].push_back(thisNCommon[iT]);
						}
	
					if (k >= 0)
					
						locMTrees[iOne][i].idProgenitor.push_back(locHalos[iOne][k].ID);
					else 
						locMTrees[iOne][i].idProgenitor.push_back(locBuffHalos[-k].ID);
						
					locMTrees[iOne][i].indexProgenitor.push_back(k);

						totCmp++;
					} else {
						totSkip++;
					}
				} // Halo Comparison
			}	// for j, k = index(j)


			/* Very important: if it turns out the halo has no likely progenitor, and has a number of particles above 
			 * minPartHalo, then we add it to the grid & the iTwo step & the iTwo grid */
			if (locMTrees[iOne][thisIndex].idProgenitor.size() == 0 && 
				locMTrees[iOne][thisIndex].mainHalo.nPart[1] > minPartHalo && iOne < iTwo)
				{
					nLocOrphans++;
					locMTrees[iOne][thisIndex].isOrphan = true;
					int addIndex = locHalos[iTwo].size();
					
					locHalos[iTwo].push_back(locMTrees[iOne][thisIndex].mainHalo);
					GlobalGrid[iTwo].AssignToGrid(locMTrees[iOne][thisIndex].mainHalo.X, addIndex);

				} else {
					locMTrees[iOne][thisIndex].isOrphan = false;
				}


		} // for i halo, the main one

		//it (locTask == 0)
			cout << "Found " << nLocOrphans << " orphan halos on task " << locTask << endl;

#endif		// ifdef ZOOM
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
	/* The trees already contain all the Halo information */
	locCleanTrees.resize(nUseCat-1);
	//allHalos.resize(nUseCat);
	//copy(locHalos[0].begin(), locHalos[0].end(), back_inserter(allHalos[iNumCat]));
};


/* This function compares the forward/backward connections to determine the unique descendant of each halo
 * TODO in full box mode, add the possibility to check the halos in the buffer and remove overlapping objects 
 * TODO Need to establish a criterion to determine to which task the halo should belong to */
void CleanTrees(int iStep)
{
	int thisIndex = 0, halosPerTask = 0, halosRemaind = 0;

        halosPerTask = int(nLocHalos[0] / totTask);
        halosRemaind = nLocHalos[0] % totTask;

	if (locTask < halosRemaind)
		halosPerTask += 1;

	if (locTask == 0)
		cout << "Cleaning Merger Tree connections, back and forth..." << endl;

	for (int kTree = 0; kTree < halosPerTask; kTree++)
	{
		int iTree = locTask + kTree * totTask;
		unsigned long long int mainID = locHalos[0][iTree].ID;

		MergerTree mergerTree;
		mergerTree.nCommon.resize(nPTypes);
		mergerTree.mainHalo = locHalos[0][iTree];
		mergerTree.isOrphan = locMTrees[0][iTree].isOrphan;
	
		/* At each step we only record the connections between halos in catalog 0 and catalog 1, without attempting at a
		 * reconstruction of the full merger history. This will be done later. */
		for (int iProg = 0; iProg < locMTrees[0][iTree].idProgenitor.size(); iProg++)
		{
			int jTree = locMTrees[0][iTree].indexProgenitor[iProg];
			unsigned long long int progID = locMTrees[0][iTree].idProgenitor[iProg];
			unsigned long long int descID = locMTrees[1][jTree].idProgenitor[0];	// The progenitor in the backwards tree
			Halo subHalo;	subHalo = locMTrees[1][jTree].mainHalo;

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
				if (locMTrees[0][iTree].isOrphan)
				{
					mergerTree.idProgenitor.push_back(locMTrees[0][iTree].mainHalo.ID);
					mergerTree.indexProgenitor.push_back(jTree);
					mergerTree.subHalos.push_back(locHalos[0][iTree]);

					for(int iT = 0; iT < nPTypes; iT++)
						mergerTree.nCommon[iT].push_back(locHalos[0][iTree].nPart[iT]);

				} else {
					mergerTree.idProgenitor.push_back(progID);
					mergerTree.indexProgenitor.push_back(jTree);
					mergerTree.subHalos.push_back(subHalo);

					for(int iT = 0; iT < nPTypes; iT++)
						mergerTree.nCommon[iT].push_back(locMTrees[0][iTree].nCommon[iT][iProg]);
				}
			}
		}	// kTree & iTree loop

		if (mergerTree.idProgenitor.size() > 0)
			locCleanTrees[iStep-1].push_back(mergerTree);

		mergerTree.Clean();
	}
	
#ifdef ZOOM
	locTreeIndex.clear();
	locTreeIndex.shrink_to_fit();
#endif
};



/* These two functions are used in mode = 1, when the MTrees are being read in from the .mtree files */
void AssignDescendant()
{
	unsigned long long int mainID = 0;
	int mainIndex = 0;

	orphanHaloIndex.clear();
	orphanHaloIndex.shrink_to_fit();

	//cout << "AssignDescendant(): " << locCleanTrees[iNumCat-1].size() << endl;

	for (int iC = 0; iC < locCleanTrees[iNumCat-1].size(); iC++)
	{
		mainID = locCleanTrees[iNumCat-1][iC].mainHalo.ID;
		mainIndex = id2Index[mainID];

		if (locCleanTrees[iNumCat-1][iC].isOrphan)
		{
			//cout << iC << " " << locCleanTrees[iNumCat-1][iC].mainHalo.ID 
				//<< " has a token progenitor: " << mainIndex << endl;
			orphanHaloIndex.push_back(mainIndex);
		} else {
			//if (locTask == 0 && iNumCat > 1 && mainIndex == 0)
			//if (iNumCat > 1 && mainIndex == 0)
			//	cout << "->" << mainID << " " << iC << " " << mainIndex << " " << locHalos[iUseCat][mainIndex].ID << endl;

			locCleanTrees[iNumCat-1][iC].mainHalo = locHalos[iUseCat][mainIndex];
		}
	}

	/* Halos have been assigned, so we can clear the map */
	id2Index.clear();	
};


void AssignProgenitor()
{
	unsigned long long int subID = 0;
	int subIndex = 0;

	orphanHaloIndex.clear();
	orphanHaloIndex.shrink_to_fit();

	//cout << "AssignDescendant(): " << locCleanTrees[iNumCat-1].size() << endl;

	for (int iC = 0; iC < locCleanTrees[iNumCat-1].size(); iC++)
	{
		//cout << iC << " " << locCleanTrees[iNumCat-1][iC].subHalos.size() << endl;

		for (int iS = 0; iS < locCleanTrees[iNumCat-1][iC].subHalos.size(); iS++ )
		{
			subID = locCleanTrees[iNumCat-1][iC].subHalos[iC].ID;
			subIndex = id2Index[subID];
			locCleanTrees[iNumCat-1][iC].subHalos[iS] = locHalos[iUseCat][subIndex];

			//if (subIndex == 0)
			//	cout << iS << ", SubIndex: " << subID << " " << subIndex << " " << endl;

		}
	}

	/* Halos have been assigned, so we can clear the map */
	//id2Index.clear();	

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



