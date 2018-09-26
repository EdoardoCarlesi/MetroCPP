/* This class takes care of input and output files, where they are located, how they are distributed across tasks and so on.
 * It uses absolute paths to input files and relies on some bash scripts to locate effectively the input catalogs and particle lists. */

#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <array>
#include <cstdio>
#include <memory>
#include <stdexcept>
//#include <stdio.h>
//#include <stdlib.h>

#include "IOSettings.h"
#include "utils.h"
#include "global_vars.h"

#define FMHIRESFAC 0.90

using namespace std;


IOSettings::IOSettings() 
{
};


IOSettings::~IOSettings() 
{
};


void IOSettings::FindCatID()
{	
	string outputSh;
	string inputSh;
	string optionsSh;
	string outputTmp;
	string lineIn;
	string cleanTmp;
	char *tmpLine;
	int iS = 0, sysOut = 0;

#ifdef ZOOM
	string boolZoom = "true";
#else
	string boolZoom = "false";
#endif

	optionsSh = pathInput + " " + haloSuffix + " " + boolZoom;

	outputTmp = pathMetroCpp + tmpIdOut;
	inputSh = pathMetroCpp + findIDsh + " " + optionsSh + " > " + outputTmp;
	cout << inputSh << endl;	

	if(FILE *f = fopen(outputTmp.c_str(), "r"))
	{
		cout << "File " << outputTmp << " found. " << endl;

	} else {

		// Execute the bash script to find out the snapshot IDs. These are usually just the snapshot numbers, 
		// but might change sometimes if there is a "hole" in between.
		cout << inputSh << endl;
		sysOut = system(inputSh.c_str());
	}

	ifstream fileIn(outputTmp);
	
	// Read the catalog IDs generated by the catalog	
	while (getline(fileIn, lineIn))
	{
		strSnaps[nCat - iS -1] = lineIn.c_str(); 
		numSnaps[nCat - iS -1] = stoi(lineIn.c_str()); 
		iS++;
	}

#ifdef CLEAN_TMP
	cleanTmp = "rm " + outputTmp;
	sysOut = system(cleanTmp.c_str());
#endif
};


void IOSettings::ReadConfigFile(string configFile)
{
	int iLine = 0;
	string lineIn, delim;
	vector<string> args; 
	
	delim = "=";

	ifstream fileCfg(configFile);

		if (!fileCfg.good())
		{
			cout << "File: " << configFile << " not found. Exiting..." << endl;
			MPI_Finalize();
			exit(0);

		} else {
#ifdef VERBOSE
			if (locTask == 0)
				cout << "Reading config file: " << configFile << endl;
#endif
		}

		while (getline(fileCfg, lineIn))
		{
			const char *lineRead = lineIn.c_str();		
	
			if (lineRead[0] != '#' && lineRead[0] != ' ')
			{
				if (locTask == 0)
				{
					args = SplitString(lineRead, delim);
	
#ifdef VERBOSE
						if (args.size() > 1) 
						{
							cout << args[0] << endl;
							cout << args[1] << endl;
						}
#endif
				}
			}

			/* The line is read in, now setup */
			if (args.size() > 1) 
				InitFromCfgFile(args);

			iLine++;

		} /* while(line) */
};


void IOSettings::InitFromCfgFile(vector<string> arg)
{
	if (arg[0] == "boxSize") 		boxSize = stof(arg[1]);
	else if (arg[0] == "haloSuffix") 	haloSuffix = arg[1];
	else if (arg[0] == "partSuffix") 	partSuffix = arg[1];
	else if (arg[0] == "haloPrefix") 	haloPrefix = arg[1];
	else if (arg[0] == "partPrefix") 	partPrefix = arg[1];
	else if (arg[0] == "pathMetroCpp") 	pathMetroCpp = arg[1];
	else if (arg[0] == "pathInput") 	pathInput = arg[1];
	else if (arg[0] == "inputFormat") 	inputFormat = arg[1];
	else if (arg[0] == "nCat")		nCat = stoi(arg[1]);
	else if (arg[0] == "nChunks")		nChunks = stoi(arg[1]);
	else if (arg[0] == "nGrid")		nGrid = stoi(arg[1]);
	else if (arg[0] == "dMaxFactor")	dMaxFactor = stof(arg[1]);
	else {
		cout << "Arg= " << arg[0] << " does not exist." << endl;
	}

	//if (arg.size() > 0)
	//	cout << arg[0] << " = " << arg[1] << endl;

	/* Just issue a warning here, in case some parameter has not been set correctly. */
	if (arg[1] == "" && locTask == 0)
		cout << "WARNING " << arg[0] << " has not been set correctly in the config file." << endl;
 
}


void IOSettings::CheckStatus()
{
	cout << "On task="	 << locTask << endl;
	cout << "nCat        = " << nCat << endl;
	cout << "nChunks     = " << nChunks << endl;
	cout << "inputFormat = " << inputFormat << endl;
	cout << "nGrid       = " << nGrid << endl;
	cout << "boxSize     = " << boxSize << endl;
	cout << "pathMetro   = " << pathMetroCpp << endl;
	cout << "pathInput   = " << pathInput << endl;
}


void IOSettings::FindCatZ()
{	
	int iZ = 0, sysOut = 0;
	string outputSh;
	string inputSh;
	string optionsSh;
	string outputTmp;
	string cleanTmp;
	string lineIn;

#ifdef ZOOM
	string boolZoom = "true";
#else
	string boolZoom = "false";
#endif

	optionsSh = pathInput + " " + haloSuffix + " " + boolZoom;
	outputTmp = pathMetroCpp + tmpZOut;
	inputSh = pathMetroCpp + findZsh + " " + optionsSh + " > " + outputTmp;
	cout << inputSh << endl;	

	if(FILE *f = fopen(outputTmp.c_str(), "r"))
	{
		cout << "File " << outputTmp << " found. " << endl;

	} else {
	
		// Execute the bash script and find out the redshifts of the snapshot files.
		// TODO this assumes AHF format! Other formats might not dump the z value in the output file
		cout << inputSh << endl;
		sysOut = system(inputSh.c_str());
	}

	ifstream fileIn(outputTmp);
		
	// Read the redshifts generated by the catalog
	while (getline(fileIn, lineIn))
	{
		const char *lineRead = lineIn.c_str();
		float thisZ = 0.0, thisA = 0.0;

		sscanf(lineRead, "%f", &thisZ);
		thisA = 1.0 / (1.0 + thisZ) ;

		redShift[nCat - iZ - 1] = thisZ;
		aFactors[nCat - iZ - 1] = thisA;

		//cout << nCat-iZ-1 << ", z = " << redShift[nCat - iZ -1] << ", a=" << aFactors[nCat - iZ -1]<< endl;
		iZ++;
	}
	
#ifdef CLEAN_TMP
	// Remove temporary files
	cleanTmp = "rm " + outputTmp;
	sysOut = system(cleanTmp.c_str());
#endif
};


void IOSettings::FindCatN()
{	
	int sysOut = 0;
	string outputSh;
	string inputSh;
	string optionsSh;
	string outputTmp;
	string lineIn;
	string cleanTmp;

#ifdef ZOOM
	string boolZoom = "true";
#else
	string boolZoom = "false";
#endif

	optionsSh = pathInput + " " + haloSuffix + " " + boolZoom;

	outputTmp = pathMetroCpp + tmpNOut;
	inputSh = pathMetroCpp + findNsh + " " + optionsSh + " > " + outputTmp;
	cout << inputSh << endl;

	if(FILE *f = fopen(outputTmp.c_str(), "r"))
	{
		cout << "File " << outputTmp << " found. " << endl;
	} else {

		sysOut = system(inputSh.c_str());
	}

		ifstream fileIn(outputTmp);
		
	// Read the number of catalogs generated by the script
	while (getline(fileIn, lineIn))
	{
		const char *lineRead = lineIn.c_str();
		sscanf(lineRead, "%d", &nCat);
	}

#ifdef CLEAN_TMP	
	// Remove temporary files
	cleanTmp = "rm " + outputTmp;
	sysOut = system(cleanTmp.c_str());
#endif
};


void IOSettings::Init()
{
	// Use only one task to read the files
	if (locTask == 0)
		FindCatN();

	// Once the catalog number has been found, communicate it to all tasks
	MPI_Bcast(&nCat, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Allocate memory for the catalog names, numbers, redshifts and a factors
	redShift.resize(nCat);
	aFactors.resize(nCat);
	strSnaps.resize(nCat);
	numSnaps.resize(nCat);

	// Now read the catalog names on the master task, broadcast everything later
	if (locTask == 0)
	{
		cout << "Found " << nCat << " redshifts in total." << endl;

		FindCatID();
		FindCatZ();
	}

	MPI_Bcast(&redShift[0], nCat, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&aFactors[0], nCat, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&numSnaps[0], nCat, MPI_INT, 0, MPI_COMM_WORLD);

	if (inputFormat == "AHF")
	{
		if (locTask == 0)
			cout << "Using " << inputFormat << " file format." << endl;

	} else if (inputFormat == "FoF") {

		if (locTask == 0)
			cout << "Format " << inputFormat << " is not supported yet. Exiting..." << endl;

		MPI_Finalize();
		exit(0);

	} else {

		if (locTask == 0)
			cout << "Format " << inputFormat << " is not supported yet. Exiting..." << endl;

		MPI_Finalize();
		exit(0);
	};

	// Convert integerts to snapshot strings on each task, it is easier than MPI_Bcast all those chars
	if (locTask != 0)
	{
		// The snapshot format is _XXX so we assume 4 char. BEWARE! Other formats might not be compatible
		char strBuff[4];

		for (int i = 0; i < nCat; i++)
		{
			sprintf(strBuff, "%03d", numSnaps[i]); 
			strSnaps[i] = strBuff;
		}
	}

};


void IOSettings::DistributeFilesAmongTasks(void)
{
	char charCpu[5], charZ[8];
	int locChunk = 0;
	string strZ;

#ifdef ZOOM
	if (locTask == 0)
		cout << "Reading halo/particle files on Task=0. Total number of tasks= " << totTask << endl; 
#else
	if (locTask == 0)
		cout << "Each task is reading " << nChunks << " halo/particle files. Total number of tasks= " << totTask << endl; 
#endif

	haloFiles.resize(nCat);
	partFiles.resize(nCat);

	for (int i = 0; i < nCat; i++)
	{
		sprintf(charZ, "%.3f", redShift[i]);	
		haloFiles[i].resize(nChunks);	/* in ZOOM mode there is only one chunk */
		partFiles[i].resize(nChunks);

		for (int j = 0; j < nChunks; j++)
		{
#ifdef ZOOM
			if (locTask == 0)
			{	
				haloFiles[i][0] = pathInput + haloPrefix + strSnaps[i] + ".z" + charZ + "." + haloSuffix;

				ifstream haloExists(haloFiles[i][j]);
				if (haloExists.fail())
				{
					haloFiles[i][0] = pathInput + haloPrefix + strSnaps[i] + ".0000.z" + charZ + "." + haloSuffix;

					ifstream haloExists(haloFiles[i][j]);
						if (haloExists.fail())
							cout << "WARNING: on task =" << locTask << " AHF_halos found as " << haloFiles[i][j] << endl;
				}

				partFiles[i][0] = pathInput + haloPrefix + strSnaps[i] + ".z" + charZ + "." + partSuffix;

				ifstream partExists(partFiles[i][j]);
				if (partExists.fail())
				{
					partFiles[i][0] = pathInput + partPrefix + strSnaps[i] + ".0000.z" + charZ + "." + partSuffix;

					ifstream partExists(partFiles[i][j]);
						if (partExists.fail())
							cout << "WARNING: on task =" << locTask << " AHF_halos found as " << partFiles[i][j] << endl;
				}

			}

#else		/* No ZOOM, distribute the files as usually */

			locChunk = j + locTask * nChunks;
			sprintf(charCpu, "%04d", locChunk);	
			haloFiles[i][j] = pathInput + haloPrefix + strSnaps[i] + "." + charCpu + ".z" + charZ + "." + haloSuffix;
			partFiles[i][j] = pathInput + haloPrefix + strSnaps[i] + "." + charCpu + ".z" + charZ + "." + partSuffix;

			ifstream haloExists(haloFiles[i][j]);
			if (haloExists.fail())
				cout << "WARNING: on task =" << locTask << " " << haloFiles[i][j] << " not found." << endl;

			ifstream partExists(partFiles[i][j]);
			if (partExists.fail())
				cout << "WARNING: on task =" << locTask << " " << haloFiles[i][j] << " not found." << endl;


#endif
			//cout << locTask << ") " << haloFiles[i][j] << endl; 
		}
	}
};


// Particle sizes have already been allocated in the ReadHalos() routines, do a safety check for the size
// TODO use read(buffer,size) to read quickly blocks of particles all at the same time 
void IOSettings::ReadParticles(void)
{
	vector<vector<unsigned long long int>> tmpParts;
	string tmpStrUrlPart, lineIn;
	const char *tmpUrlPart;
	unsigned long long int locHaloID = 0, partID = 0;
	unsigned int iTmpParts = 0, iLocParts = 0, iLine = 0, nPartHalo = 0;
	unsigned int nFileHalos = 0, iLocHalos = 0, iTmpHalos = 0;
	int partType = 0;

	tmpParts.resize(nPTypes);
	locParts[iUseCat].resize(nLocHalos[iUseCat]);

#ifdef VERBOSE
	cout << locTask << ") Reading particles for n halos = " << nLocHalos[iUseCat] << endl;
#endif

#ifdef ZOOM	/* Only read on one task */
	if (locTask == 0)
	{
#endif

	for (int iChunk = 0; iChunk < nChunks; iChunk++)
	{
		tmpUrlPart = partFiles[iNumCat][iChunk].c_str();
		ifstream fileIn(tmpUrlPart);

		/* Reset temporary variables */
		iTmpParts = 0;
		iTmpHalos = 0;
		iLine = 0;

		if (!fileIn.good())
		{
			cout << "File: " << tmpUrlPart << " not found on task=" << locTask << endl;
		} else {
			if (locTask == 0)
	        		cout << "Reading particle file: " << tmpUrlPart << endl;
		}

		while (getline(fileIn, lineIn))
		{
			const char *lineRead = lineIn.c_str();		
			
			if (iLine == 0)
			{
				if (inputFormat == "AHF")
		         	       sscanf(lineRead, "%d", &nFileHalos);

				iLine++;

			} else if (iLine == 1) {

				if (inputFormat == "AHF")
		        	        sscanf(lineRead, "%u %llu", &nPartHalo, &locHaloID);
					//cout << locTask << " " << nPartHalo << " " << iLine << " " << locHaloID << endl;
#ifdef ZOOM
			if (locHalos[iUseCat][iLocHalos].ID == locHaloID)
#endif
				locParts[iUseCat][iLocHalos].resize(nPTypes);

				iLine++;
			} else {

				if (inputFormat == "AHF")
		        	        sscanf(lineRead, "%llu %d", &partID, &partType);

				//cout << locTask << " " << iLine << " " << partID << endl;
				tmpParts[partType].push_back(partID);
				iTmpParts++;
				iLocParts++;

#ifdef ZOOM
				if (iTmpParts == nPartHalo)
#else
				if (iTmpParts == locHalos[iUseCat][iLocHalos].nPart[nPTypes])
#endif
				{	

					// Sort the ordered IDs
					for (int iT = 0; iT < nPTypes; iT++)
					{	
						if (tmpParts[iT].size() > 0)
						{
#ifdef ZOOM
							if (locHalos[iUseCat][iLocHalos].ID == locHaloID)
							{
#endif
								//cout << locTask << " " << tmpParts[iT].size() << endl;
								sort(tmpParts[iT].begin(), tmpParts[iT].end());
						
								locParts[iUseCat][iLocHalos][iT].insert(
										locParts[iUseCat][iLocHalos][iT].end(), 
											tmpParts[iT].begin(), tmpParts[iT].end());

								//cout << locTask << " " << locParts[iUseCat][iLocHalos][iT][0] << endl;
#ifdef ZOOM
							} // Zoom mode, making sure the current halo is in the list of the high-res ones
#endif
							// Clean the temporary read-in buffer
							tmpParts[iT].clear();
							tmpParts[iT].shrink_to_fit();
						}
					}

#ifdef ZOOM
					if (locHalos[iUseCat][iLocHalos].ID == locHaloID)
#endif
						iLocHalos++;

					iLocParts++;

					// Set/reset some counters
					iTmpParts = 0;
					iTmpHalos++;

					// Check if all of the halos in the current file chunk have been read in
					if (iTmpHalos == nFileHalos)
					{
						iTmpHalos = 0;
						iLine = 0;	
					} else {
						iLine = 1;	
					}
				}
			} // else iLine not 0 or 1
		} // end while
	} // End for loop on file chunks

#ifdef ZOOM
	/* There is no actual loop but we close the reading of the file on task 0 */

#ifdef VERBOSE
	} else {
		cout << "Task=" << locTask << " is waiting for communication from Task 0 " << endl;
#endif
	}	
#endif

	
#ifdef VERBOSE
	cout << "All particle files for " << iLocHalos << " halos have been read read on task " << locTask << endl;
#endif
};
 

/* Using AHF by default */
void IOSettings::ReadHalos()
{
	unsigned int nPartHalo = 0, nPTypes = 0, iTmpHalos = 0, nTmpHalos = 0, iChunk = 0; 
	const char *tmpUrlHalo, *lineHead = "#";
	string tmpStrUrlHalo, lineIn;
	vector<Halo> tmpHalos;

	// Initialize outside of the files loop
	int iLocHalos = 0; 

#ifdef ZOOM	/* Only read on one task */
	if (locTask == 0)
	{
#endif

	for (int iChunk = 0; iChunk < nChunks; iChunk++)
	{
		tmpUrlHalo = haloFiles[iNumCat][iChunk].c_str();
		nTmpHalos = NumLines(tmpUrlHalo);	
		tmpHalos.resize(nTmpHalos); 
		iTmpHalos = 0;

		ifstream fileIn(tmpUrlHalo);
	
		if (!fileIn.good())
		{
			cout << "File: " << tmpUrlHalo << " not found on task=" << locTask << endl;
		} else { 
			if (locTask == 0)
	       			cout << "Reading " << nTmpHalos << " halos from file: " << tmpUrlHalo << endl;
		}

		while (getline(fileIn, lineIn))
		{
			const char *lineRead = lineIn.c_str();		
		
			if (lineRead[0] != lineHead[0])
			{
				ReadLineAHF(lineRead, &tmpHalos[iTmpHalos]);
				nPartHalo = tmpHalos[iTmpHalos].nPart[nPTypes];	// All particle types!

#ifndef ZOOM
				// Assign halo to its nearest grid point - assign the absolute local index number
				// Halos on the local chunk have POSITIVE index, halos on the buffer NEGATIVE 
				GlobalGrid[iUseCat].AssignToGrid(tmpHalos[iTmpHalos].X, iLocHalos);
				locId2Index.insert(pair <unsigned long long int, int> (tmpHalos[iTmpHalos].ID, iLocHalos));
#else
				id2Index.insert(pair <unsigned long long int, int> (tmpHalos[iTmpHalos].ID, iLocHalos));
#endif
				iLocHalos++;
				iTmpHalos++;
			}
		}	/* While Read Line */

		fileIn.close();

#ifdef ZOOM	/* When in ZOOM mode ONLY! store the high density region halos */
		iLocHalos = 0;

		for (int iH = 0; iH < tmpHalos.size(); iH++)
		{
			if (tmpHalos[iH].fMhires >= FMHIRESFAC)
			{
				locHalos[iUseCat].push_back(tmpHalos[iH]);
				iLocHalos++;
			}
		}

		//if (locTask == 0)
		//	cout << "Total halos: " << tmpHalos.size() << ", high-res= " << iLocHalos << endl;

		tmpHalos.clear();
		tmpHalos.shrink_to_fit();
#else
	
#ifdef VERBOSE
		cout << "NHalos: " << tmpHalos.size() << " on task=" << locTask << endl;
#endif
		/* Append to the locHalo file */
		locHalos[iUseCat].insert(locHalos[iUseCat].end(), tmpHalos.begin(), tmpHalos.end());
		tmpHalos.clear();
		tmpHalos.shrink_to_fit();
#endif
	} // Loop on files per task
	
#ifdef ZOOM
	}
#endif

	nLocHalos[iUseCat] = iLocHalos;

#ifndef ZOOM
	// After reading in all the catalogs, find out, sort and remove duplicates of nodes being allocated to the task
	GlobalGrid[iUseCat].SortLocNodes();
#endif

#ifdef VERBOSE
	cout << "On task=" << locTask << ", " << nLocHalos[iUseCat] * sizeHalo /1024/1024 << " MB haloes read " << endl; 
#endif
};


void IOSettings::ReadLineAHF(const char * lineRead, Halo *halo)
{
	float dummy, vHalo;
	unsigned int tmpNpart, nGas = 0, nStar = 0;

	/* AHF file structure:
	   ID(1)  hostHalo(2)     numSubStruct(3) Mvir(4) npart(5)        Xc(6)   Yc(7)   Zc(8)   VXc(9)  VYc(10) VZc(11) 
	   Rvir(12)        Rmax(13)        r2(14)  mbp_offset(15)  com_offset(16)  Vmax(17)        v_esc(18)       sigV(19)
           lambda(20)      lambdaE(21)     Lx(22)  Ly(23)  Lz(24)  b(25)   c(26)   Eax(27) Eay(28) Eaz(29) Ebx(30) Eby(31) 
	   Ebz(32) Ecx(33) Ecy(34) Ecz(35) ovdens(36)      nbins(37)       fMhires(38)     Ekin(39)        Epot(40)        
	   SurfP(41)       Phi0(42)        cNFW(43)        n_gas(44)       M_gas(45)       lambda_gas(46)  lambdaE_gas(47) 
	   Lx_gas(48)      Ly_gas(49)      Lz_gas(50)      b_gas(51)       c_gas(52)       Eax_gas(53)     Eay_gas(54)     
	   Eaz_gas(55)     Ebx_gas(56)     Eby_gas(57)     Ebz_gas(58)     Ecx_gas(59)     Ecy_gas(60)     Ecz_gas(61)     
	   Ekin_gas(62)    Epot_gas(63)    n_star(64)      M_star(65)      lambda_star(66) lambdaE_star(67)        
           Lx_star(68)     Ly_star(69)     Lz_star(70)     b_star(71)      c_star(72)      Eax_star(73)    Eay_star(74) Eaz_star(75)    
	   Ebx_star(76)    Eby_star(77)    Ebz_star(78)    Ecx_star(79)    Ecy_star(80)    Ecz_star(81)    Ekin_star(82)Epot_star(83) */

	/* Col:		    1   2    3  4  5  6  7  8  9 10 11 */
	sscanf(lineRead, "%llu %llu %d %f %d \
			  %f   %f   %f %f %f %f \
			  %f   %f   %f %f %f %f %f %f %f %f \
			  %f   %f   %f \
			  %f   %f   %f %f %f %f %f %f %f %f \
			  %f   %f   %f %f %f %f %f %f %f %f ",
 
			&halo->ID, &halo->hostID, &halo->nSub, &halo->mTot, &tmpNpart, 
			&halo->X[0], &halo->X[1], &halo->X[2], &halo->V[0], &halo->V[1], &halo->V[2], 				// 11
			&halo->rVir, &dummy, &halo->rsNFW, &dummy, &dummy, &halo->vMax, &dummy, &halo->sigV, &halo->lambda, &dummy, // 21
			&halo->L[0], &halo->L[1], &halo->L[2],									// 24
			&dummy, &dummy, &dummy, &dummy,   &dummy, &dummy, &dummy, &dummy, &dummy, &dummy,			// 34
			&dummy, &dummy, &dummy, &halo->fMhires, &dummy, &dummy, &dummy, &dummy, &halo->cNFW, &dummy);		 // 44

	/* Particle numbers were not allocated correctly sometimes, so let's reset them carefully */
	nGas = 0; nStar = 0;
	halo->nPart[0] = nGas;
	halo->nPart[1] = tmpNpart - nGas - nStar;
	halo->nPart[2] = 0;
	halo->nPart[3] = nStar;
	halo->nPart[4] = 0;
	halo->nPart[5] = 0;
	halo->nPart[6] = tmpNpart;

	/* Compute max velocity and sub box edges while reading the halo file ---> this is used to compute the buffer zones */
	vHalo = VectorModule(halo->V);

	if (vHalo > locVmax)
		locVmax = vHalo;
};


   
void IOSettings::WriteTrees()
{

};

