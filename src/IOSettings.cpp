/*
 *   METROC++: MErger TRees On C++, a scalable code for the computation of merger trees in cosmological simulations.
 *   Copyright (C) Edoardo Carlesi 2018-2019
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


/* 
 * IOSettings.cpp:
 * This class takes care of input and output files, where they are located, how they are distributed across tasks and so on.
 * It uses absolute paths to input files and relies on some bash scripts to locate effectively the input catalogs and particle lists. 
 */

#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <array>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <map>

#include "Cosmology.h"
#include "IOSettings.h"
#include "utils.h"
#include "spline.h"
#include "global_vars.h"

#ifdef IDADD
#define idADD 1000000000
#endif

using namespace std;


IOSettings::IOSettings() 
{
};


IOSettings::~IOSettings() 
{
};


void IOSettings::SetCosmology(Cosmology *Cosmo)
{
	if (cosmologicalModel == "WMAP7")
	{
		if (locTask == 0)
			cout << "Using WMAP7 cosmology." << endl;
			
		pathA  = pathMetroCpp + dataAWMAP7;
		pathPk = pathMetroCpp + dataPkWMAP7;
		Cosmo->SetWMAP7();

	} else if (cosmologicalModel == "Planck") {

		if (locTask == 0)
			cout << "Using Planck cosmology." << endl;
			
		pathA  = pathMetroCpp + dataAPlanck;
		pathPk = pathMetroCpp + dataPkPlanck;
		Cosmo->SetPlanck();

	} else {

		if (locTask == 0)
		{
			cout << "Error. Cosmological model << " << cosmologicalModel << ">> is not implemented. " << endl; 
			cout << "Try setting << WMAP7 >> or << Planck >> instead." << endl;
			cout << "Exiting..." << endl;
		}

		MPI_Finalize();
		exit(0);
	}

	Cosmo->a  = ReadA();
	Cosmo->pk = ReadPk();
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
	else if (arg[0] == "outSuffix")		outSuffix = arg[1];
	else if (arg[0] == "haloPrefix") 	haloPrefix = arg[1];
	else if (arg[0] == "partPrefix") 	partPrefix = arg[1];
	else if (arg[0] == "outPrefix")		outPrefix = arg[1];
	else if (arg[0] == "cpuString")		cpuString = arg[1];
	else if (arg[0] == "splitString")	splitString = arg[1];
	else if (arg[0] == "pathMetroCpp") 	pathMetroCpp = arg[1];
	else if (arg[0] == "pathInput") 	pathInput = arg[1];
	else if (arg[0] == "inputFormat") 	inputFormat = arg[1];
	else if (arg[0] == "nSnapsUse")		nSnapsUse = stoi(arg[1]);
	else if (arg[0] == "nChunks")		nChunks = stoi(arg[1]);
	else if (arg[0] == "nGrid")		nGrid = stoi(arg[1]);
	else if (arg[0] == "facOrphanSteps")	facOrphanSteps = stoi(arg[1]);
	else if (arg[0] == "maxOrphanSteps")	maxOrphanSteps = stoi(arg[1]);
	else if (arg[0] == "minPartHalo")	minPartHalo = stoi(arg[1]);
	else if (arg[0] == "minPartCmp")	minPartCmp = stoi(arg[1]);
	else if (arg[0] == "pathOutput")	pathOutput = arg[1];
	else if (arg[0] == "nTreeChunks")	nTreeChunks = stoi(arg[1]);
	else if (arg[0] == "cosmologicalModel")	cosmologicalModel = arg[1];
	else cout << "Arg= " << arg[0] << " is useless or redundant and will be ignored." << endl;

	/* Just issue a warning here, in case some parameter has not been set correctly. */
	if (arg[1] == "" && locTask == 0)
		cout << "WARNING " << arg[0] << " has not been set correctly in the config file." << endl;
 
#ifndef ZOOM
	if (arg[0] == "nGrid" && locTask == 0)
		if (nGrid < 1)
		{
			cout << "ERROR: " << endl; 
			cout << "Cannot run with parameter 'nGrid = 0' in Full Box mode." << endl;
			cout << "Please set nGrid to a positive int value in the configuration file." << endl;
			cout << "Exiting program..." << endl;
			exit(0);
		}
#endif

	/* Sanity check: do the specified folders exist? */
	if (arg[0] == "pathInput" && locTask == 0)
		CheckPath(pathInput);

	if (arg[0] == "pathMetroCpp" && locTask == 0)
		CheckPath(pathMetroCpp);
	
	if (arg[0] == "pathOutput" && locTask == 0)
		CheckPath(pathOutput);

	if (arg[0] == "nChunks" && locTask == 0)
		if (nChunks < 1)
		{
			cout << "ERROR: " << endl; 
			cout << "Cannot run with parameter 'nChunks = 0'. Set nChunks to a positive int value in the configuration file." << endl;
			cout << "Exiting program..." << endl;
			MPI_Finalize();
			exit(0);
		}

}



void IOSettings::CheckStatus()
{
	cout << "On task     = " << locTask << endl;
	cout << "nSnaps      = " << nSnaps << endl;
	cout << "nChunks     = " << nChunks << endl;
	cout << "inputFormat = " << inputFormat << endl;
	cout << "nGrid       = " << nGrid << endl;
	cout << "boxSize     = " << boxSize << endl;
	cout << "pathMetro   = " << pathMetroCpp << endl;
	cout << "pathInput   = " << pathInput << endl;

	for (int i=0; i<5; i++)
	{
		locHalos[0][i].Info();
		locHalos[1][i].Info();
		cout << locParts[0][i].size() <<  " type " << locParts[0][i][0].size() << endl;
	}
}


void IOSettings::FindCatID()
{	
	int iS = 0, sysOut = 0;
	string outputSh;
	string inputSh;
	string optionsSh;
	string outputTmp;
	string lineIn;
	string cleanTmp;

	optionsSh = pathInput + " " + haloSuffix + " " + cpuString + " " + splitString + " " + haloPrefix;

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
	
	/* Read the catalog IDs generated by the catalog */
	while (getline(fileIn, lineIn))
	{
		strSnaps[nSnaps - iS -1] = lineIn.c_str(); 
		numSnaps[nSnaps - iS -1] = stoi(lineIn.c_str()); 
		iS++;
	}

#ifdef CLEAN_TMP
	cleanTmp = "rm " + outputTmp;
	sysOut = system(cleanTmp.c_str());
#endif
};


void IOSettings::FindCatZ()
{	
	int iZ = 0, sysOut = 0;
	string outputSh;
	string inputSh;
	string optionsSh;
	string outputTmp;
	string cleanTmp;
	string lineIn;

	optionsSh = pathInput + " " + haloSuffix + " " + cpuString + " " + haloPrefix;
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

		redShift[nSnaps - iZ - 1] = thisZ;
		aFactors[nSnaps - iZ - 1] = thisA;

		//cout << nSnaps-iZ-1 << ", z = " << redShift[nSnaps - iZ -1] << ", a=" << aFactors[nSnaps - iZ -1]<< endl;
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

	optionsSh = pathInput + " " + haloSuffix + " " + cpuString + " " + splitString;

	outputTmp = pathMetroCpp + tmpNOut;
	inputSh = pathMetroCpp + findNsh + " " + optionsSh + " > " + outputTmp;
	cout << inputSh << endl;

	if(FILE *f = fopen(outputTmp.c_str(), "r"))
	{
		cout << "File " << outputTmp << " found. " << endl;
	} else {
		//if (locTask == 0)
		//	cout << "Executing script: " << inputSh << endl;

		sysOut = system(inputSh.c_str());
	}

	ifstream fileIn(outputTmp);
	
	//	WARNING TODO this number is set manually in the .cfg file!!!!	
	// Read the number of catalogs generated by the script
	while (getline(fileIn, lineIn))
	{
		const char *lineRead = lineIn.c_str();
		sscanf(lineRead, "%d", &nSnaps);
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
		FindCatN();	// FIXME the number of catalogs is read in from the cfg file

	// Once the catalog number has been found, communicate it to all tasks
	MPI_Bcast(&nSnaps, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Allocate memory for the catalog names, numbers, redshifts and a factors
	redShift.resize(nSnaps);
	aFactors.resize(nSnaps);
	strSnaps.resize(nSnaps);
	numSnaps.resize(nSnaps);

	// Now read the catalog names on the master task, broadcast everything later
	if (locTask == 0)
	{
		cout << "Found " << nSnaps << " redshifts in total." << endl;

		FindCatID();
		FindCatZ();
	}

	MPI_Bcast(&redShift[0], nSnaps, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&aFactors[0], nSnaps, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&numSnaps[0], nSnaps, MPI_INT, 0, MPI_COMM_WORLD);

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

		for (int iF = 0; iF < nSnaps; iF++)
		{
			sprintf(strBuff, "%03d", numSnaps[iF]); 
			strSnaps[iF] = strBuff;
		}
	}

};



void IOSettings::DistributeFilesAmongTasks(void)
{
	int nLocRemind = 0, locChunk = 0;
	char charCpu[5], charZ[8];
	string strZ;

#ifdef ZOOM
	nLocChunks = nChunks;

	if (locTask == 0)
		cout << "Reading halo and particle files on Task=0. Total number of tasks= " << totTask << endl; 
#else
	nLocChunks = int (nChunks / totTask);
	nLocRemind = nChunks % totTask;		//FIXME check this setting
#endif

	haloFiles.resize(nSnaps);
	partFiles.resize(nSnaps);

	for (int iF = 0; iF < nSnaps; iF++)
	{
		sprintf(charZ, "%.3f", redShift[iF]);	

#ifdef ZOOM
		haloFiles[iF].resize(nLocChunks);	/* in ZOOM mode there is (usually) only one chunk */
		partFiles[iF].resize(nLocChunks);
#endif

		for (int jF = 0; jF < nLocChunks; jF++)
		{
#ifdef ZOOM
			if (locTask == 0)
			{	
				haloFiles[iF][0] = pathInput + haloPrefix + strSnaps[iF] + ".z" + charZ + "." + haloSuffix;
				//cout << "TEST " << haloFiles[iF][0] << " chunks: " << nLocChunks << endl;

				ifstream haloExists(haloFiles[iF][jF]);
				if (haloExists.fail())
				{
					haloFiles[iF][0] = pathInput + haloPrefix + strSnaps[iF] + ".0000.z" + charZ + "." + haloSuffix;

					ifstream haloExists(haloFiles[iF][jF]);
						if (haloExists.fail())
						{
							cout << "ERROR on task =" << locTask 
								<< " AHF_halos file not found. " << haloFiles[iF][0] << endl;
							exit(0);		
						}
				}

				partFiles[iF][0] = pathInput + haloPrefix + strSnaps[iF] + ".z" + charZ + "." + partSuffix;
				//cout << "---------> TEST " << partFiles[iF][0] << endl;

				ifstream partExists(partFiles[iF][jF]);
				if (partExists.fail())
				{
					partFiles[iF][0] = pathInput + partPrefix + strSnaps[iF] + ".0000.z" + charZ + "." + partSuffix;

					ifstream partExists(partFiles[iF][jF]);
						if (partExists.fail())
							cout << "WARNING: on task =" << locTask << " AHF_halos found as " << partFiles[iF][jF] << endl;
				}

			}

#else			/* No ZOOM, distribute the files as usually */

			locChunk = jF + locTask * nLocChunks;

			if (locChunk < nChunks) // Make sure we are not looking for catalog chunks beyond the boundaries
			{

				sprintf(charCpu, "%04d", locChunk);	
				string locHaloFile = pathInput + haloPrefix + strSnaps[iF] + "." + charCpu + ".z" + charZ + "." + haloSuffix;
				string locPartFile = pathInput + haloPrefix + strSnaps[iF] + "." + charCpu + ".z" + charZ + "." + partSuffix;

				//cout << "On Task= " << locTask << " Halo:" << locHaloFile << " Part:" << locPartFile << endl;
	
				haloFiles[iF].push_back(locHaloFile);
				partFiles[iF].push_back(locPartFile);

				ifstream haloExists(haloFiles[iF][jF]);
				if (haloExists.fail())
				{
					cout << "ERROR: on task =" << locTask << " " << haloFiles[iF][jF] << " not found." << endl;
					exit(0);
				}

				ifstream partExists(partFiles[iF][jF]);
				if (partExists.fail())
				{
					cout << "WARNING: on task =" << locTask << " " << haloFiles[iF][jF] << " not found." << endl;
					exit(0);
				}
			}

#endif
			//cout << locTask << ") " << haloFiles[i][j] << endl; 
		}
	}
};


/* Particle sizes have already been allocated in the ReadHalos() routines, do a safety check for the size */
void IOSettings::ReadParticles(void)
{
	int iTmpParts = 0, iLocParts = 0, iLine = 0, nPartHalo = 0, partType = 0, iPartMulti = 0;
	unsigned int nFileHalos = 0, iLocHalos = 0, iTmpHalos = 0;
	uint64_t locHaloID = 0, partID = 0;
	vector<vector<uint64_t>> tmpParts;
	string tmpStrUrlPart, lineIn;
	const char *tmpUrlPart;

#ifdef VERBOSE
	cout << "onTask=" << locTask << " part size: " << locParts[iUseCat].size() << endl;
#endif

	tmpParts.resize(nPTypes);

	locParts[iUseCat].resize(nLocHalos[iUseCat]);
	nLocChunks = haloFiles[iNumCat].size();

	//cout << locTask << ") Reading particles for n halos = " << nLocHalos[iUseCat] << " nP: " << locParts[iUseCat].size() << endl;

#ifdef VERBOSE
	cout << locTask << ") Reading particles for n halos = " << nLocHalos[iUseCat] << endl;
#endif

#ifdef ZOOM	/* Only read on one task */
	if (locTask == 0)
	{
#endif

	for (int iChunk = 0; iChunk < nLocChunks; iChunk++)
	{
		tmpUrlPart = partFiles[iNumCat][iChunk].c_str();
		ifstream fileIn(tmpUrlPart);

		/* Reset temporary variables */
		iTmpParts = 0;
		iTmpHalos = 0;
		iLine = 0;

		if (!fileIn.good())
		{
			cout << "ERROR: File " << tmpUrlPart << " not found on task=" << locTask << endl;
			exit(0);
		} else {
			if (locTask == 0 && iChunk == 0)
	        		cout << "Reading particle file: " << tmpUrlPart << endl;
		}

		while (getline(fileIn, lineIn))
		{
			const char *lineRead = lineIn.c_str();		
			
			if (iLine == 0)
			{
		         	sscanf(lineRead, "%d", &nFileHalos);
				iLine++;

			} else if (iLine == 1) {

		        	sscanf(lineRead, "%u %lu", &nPartHalo, &locHaloID);
#ifdef IDADD
				if (locHaloID < idADD)
					locHaloID += idADD;
#endif

#ifdef ZOOM
			if (locHalos[iUseCat][iLocHalos].ID == locHaloID)
#endif
				locParts[iUseCat][iLocHalos].resize(nPTypes);
				iLine++;
			} else {
#ifdef NOPTYPE
	        	        sscanf(lineRead, "%lu", &partID);
				partType = 1;
#else	
	        	        sscanf(lineRead, "%lu %d", &partID, &partType);
#endif
				Particle thisParticle;
				thisParticle.haloID = locHaloID;
				thisParticle.type   = partType;
		
				locMapParts[iUseCat][partID].push_back(thisParticle);
		
				if (locMapParts[iUseCat][partID].size() > 1)
					iPartMulti++;

				tmpParts[partType].push_back(partID);
				iTmpParts++;
				iLocParts++;

				if (iTmpParts == nPartHalo)
				{	

					/* Sort the ordered IDs */
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

//cout << iUseCat << " N particles: " << locMapParts[iUseCat].size() << " iLocParts: " << iLocParts << " Duplicates: " << iPartMulti
//		<< " total: " << locMapParts[iUseCat].size() + iPartMulti << endl;
#ifdef VERBOSE
	cout << " N particles: " << locMapParts[iUseCat].size() << " iLocParts: " << iLocParts << " Duplicates: " << iPartMulti
		<< " total: " << locMapParts[iUseCat].size() + iPartMulti << endl;
#endif
};
 

/* Using AHF by default */
void IOSettings::ReadHalos()
{
	unsigned int nPartHalo = 0, nPTypes = 0, iTmpHalos = 0, nTmpHalos = 0, iChunk = 0, iLocHalos = 0; 
	const char *tmpUrlHalo, *lineHead = "#";
	string tmpStrUrlHalo, lineIn;
	vector<Halo> tmpHalos;

#ifdef ZOOM	/* Only read on one task */
	if (locTask == 0)
	{
#endif

#ifdef VERBOSE
	cout << "onTask=" << locTask << " halo size: " << locHalos[iUseCat].size() << endl;
#endif

	nLocChunks = haloFiles[iNumCat].size();
	//cout << locTask << ", " << iNumCat << ", " << nLocChunks << endl;

	for (int iChunk = 0; iChunk < nLocChunks; iChunk++)
	{
		tmpUrlHalo = haloFiles[iNumCat][iChunk].c_str();
		nTmpHalos = NumLines(tmpUrlHalo);	
		tmpHalos.resize(nTmpHalos); 
		iTmpHalos = 0;

		ifstream fileIn(tmpUrlHalo);
	
		if (!fileIn.good())
		{
			cout << "File: " << tmpUrlHalo << " not found on task=" << locTask << endl;
			exit(0);
		} else { 
			if (locTask == 0 && iChunk == 0)
	       			cout << "Reading " << nTmpHalos << " halos from file: " << tmpUrlHalo << endl;
		}

		while (getline(fileIn, lineIn))
		{
			const char *lineRead = lineIn.c_str();		
		
			if (lineRead[0] != lineHead[0])
			{
#ifdef AHF_CB
				ReadLineAHF_CB(lineRead, &tmpHalos[iTmpHalos]);
#else
				ReadLineAHF(lineRead, &tmpHalos[iTmpHalos]);
#endif

				nPartHalo = tmpHalos[iTmpHalos].nPart[nPTypes];	// All particle types!
#ifndef ZOOM
				// Assign halo to its nearest grid point - assign the absolute local index number
				// Halos on the local chunk have POSITIVE index, halos on the buffer NEGATIVE 
				GlobalGrid[iUseCat].AssignToGrid(tmpHalos[iTmpHalos].X, iLocHalos);
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
			locHalos[iUseCat].push_back(tmpHalos[iH]);
			iLocHalos++;
		}

		tmpHalos.clear();
		tmpHalos.shrink_to_fit();
#else
	
#ifdef VERBOSE
		cout << "NHalos: " << tmpHalos.size() << " on task=" << locTask << endl;
#endif
		/* Append to the locHalo file */
		//locHalos[iUseCat].insert(locHalos[iUseCat].end(), tmpHalos.begin(), tmpHalos.end()); // TODO maybe this works but
												       // let's stick to the safe side
		for (int iH = 0; iH < tmpHalos.size(); iH++)
			locHalos[iUseCat].push_back(tmpHalos[iH]);

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
};


/* Read a set of previously computed merger trees, these will be post processed and smoothed */
void IOSettings::ReadTrees()
{
	uint64_t hostHaloID = 0, progHaloID = 0;
	int hostPart = 0, progPart = 0, orphanHalo = 0; 
	int iLine = 0, thisNumCat = 0, nProgHalo = 0;
	int commPart = 0;

	char charChunk[2], charSnap[4];
	const char *lineHead;
	string urlTree, lineIn;

	thisNumCat = iNumCat-1; 
	sprintf(charSnap, "%03d", thisNumCat);	
	sprintf(charChunk, "%d", locTask);	

	lineHead = "#";

	if (nTreeChunks != totTask)
	{
		if (locTask == 0) 
			cout << "ERROR: nTreeChunks in the .cfg file has to be equal to the number of MPI Tasks in use. " << endl;

		MPI_Finalize();
		exit(0);
		
	} else {
		urlTree = pathTree + outPrefix + charSnap + "." + charChunk + ".mtree";
		ifstream fileIn(urlTree);

		if (!fileIn.good())
		{
			cout << "ERROR: File " << urlTree << " not found on task=" << locTask << endl;
			MPI_Finalize();
			exit(0);
		} else {
			if (locTask == 0)
	       			cout << "Reading tree file: " << urlTree << endl;
		}	

		MergerTree mergerTree;

		while (getline(fileIn, lineIn))
		{
			const char *lineRead = lineIn.c_str();

			if (lineRead[0] != lineHead[0])
			{
				if (iLine == 0)
				{
		        		sscanf(lineRead, "%lu  %d  %d  %d", &hostHaloID, &hostPart, &nProgHalo, &orphanHalo);

					mergerTree.mainHalo.ID = hostHaloID;
					mergerTree.mainHalo.nPart[1] = hostPart;	//TODO this assumes n tot particles = n DM

					mergerTree.nCommon.resize(nPTypes);
					
					for (int iC = 0; iC < nPTypes; iC++)
						mergerTree.nCommon[iC].resize(nProgHalo);

					mergerTree.progHalo.resize(nProgHalo);
					mergerTree.idProgenitor.resize(nProgHalo);

					if (orphanHalo == 1)
					{
						mergerTree.isOrphan = true;
						mergerTree.progHalo[0].isToken = true;
					} else {
						mergerTree.isOrphan = false;
						mergerTree.progHalo[0].isToken = false;
					}

					iLine++;
				} 
				else if (iLine > 0 && iLine < nProgHalo+1)	/* Read-in properties of progenitors */
				{
		        		sscanf(lineRead, "%d  %lu  %d", &commPart, &progHaloID, &progPart);

					mergerTree.idProgenitor[iLine-1] = progHaloID;
					mergerTree.nCommon[1][iLine-1] = commPart;
					mergerTree.progHalo[iLine-1].ID = progHaloID;
					mergerTree.progHalo[iLine-1].nPart[1] = progPart;
					iLine++;
				}

				/* Finished reading in progenitor halos */
				if (iLine == nProgHalo+1)
				{
					locCleanTrees[iNumCat-1].push_back(mergerTree);
					mergerTree.Clean();
					iLine = 0;
				}
			}
		}
	}	/* End if using the correct number of tasks */

};



void IOSettings::ReadLineAHF(const char * lineRead, Halo *halo)
{
	float dummy, vHalo;
	unsigned int tmpNpart = 0, nGas = 0, nStar = 0;
	//uint64_t dummyID;
	int dummyID;

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
	sscanf(lineRead, "%lu  %d   %d %f %d \
			  %f   %f   %f %f %f %f \
			  %f   %f   %f %f %f %f %f %f %f %f \
			  %f   %f   %f \
			  %f   %f   %f %f %f %f %f %f %f %f \
			  %f   %f   %f %f %f %f %f %f %f %f ",
 
			&halo->ID,   &dummyID, &halo->nSub, &halo->mTot, &tmpNpart, 
			&halo->X[0], &halo->X[1], &halo->X[2], &halo->V[0], &halo->V[1], &halo->V[2], 				// 11
			&halo->rVir, &dummy, &halo->rsNFW, &dummy, &dummy, &halo->vMax, &dummy, &halo->sigV, &halo->lambda, &dummy, // 21
			&halo->L[0], &halo->L[1], &halo->L[2],									// 24
			&dummy, &dummy, &dummy, &dummy,   &dummy, &dummy, &dummy, &dummy, &dummy, &dummy,			// 34
			&dummy, &dummy, &dummy, &halo->fMhires, &dummy, &dummy, &dummy, &dummy, &halo->cNFW, &dummy);		 // 44

	/* Particle numbers were not allocated correctly sometimes, so let's reset them carefully */
	nGas = 0; nStar = 0;

/*
#ifdef IDADD
	if (dummyID < idADD)
		dummyID = dummyID + idADD;
#endif

	halo->ID = dummyID;
*/

#ifdef NOPTYPE
	halo->nPart[0] = 0;
	halo->nPart[1] = tmpNpart;
	
	/*if (tmpNpart > 10000)
	{ cout << lineRead << endl;
		cout << halo->ID << " " << tmpNpart << endl;}*/
#else
	halo->nPart[0] = nGas;
	halo->nPart[1] = tmpNpart - nGas - nStar;

	if (nPTypes > 2)
	{
		halo->nPart[2] = 0;
		
		if (nPTypes > 3)
		{
			halo->nPart[3] = 0;
		
			if (nPTypes > 4)
			{
				halo->nPart[4] = nStar;
					
					if (nPTypes > 5)
					{
						halo->nPart[5] = 0;
					}
			}
		}
	}
#endif
	/* Compute max velocity and sub box edges while reading the halo file ---> this is used to compute the buffer zones */
	vHalo = VectorModule(halo->V);

	if (vHalo > locVmax)
		locVmax = vHalo;
};

#ifdef AHF_CB
// TODO!!!!!! Enable a different output format --> CB format does not have halo IDs, which is problematic for the MergerTree algorithm
void IOSettings::ReadLineAHF_CB(const char * lineRead, Halo *halo, int lineID)
{
	float dummy, vHalo, dummyID;
	unsigned int tmpNpart, nGas = 0, nStar = 0;

	/* AHF file structure:
	 npart(1)       fMhires(2)      Xc(3)   Yc(4)   Zc(5)   VXc(6)  VYc(7)  VZc(8)  Mvir(9) Rvir(10)        Vmax(11)        Rmax(12)
  	 sigV(13)        lambda(14)      Lx(15)  Ly(16)  Lz(17)  a(18)   Eax(19) Eay(20) Eaz(21) b(22)   
	 Ebx(23) Eby(24) Ebz(25) c(26)   Ecx(27) Ecy(28) Ecz(29) ovdens(30)      Redge(31)       
	 nbins(32)       Ekin(33)        Epot(34)        mbp_offset(35)  com_offset(36)  r2(37)  lambdaE(38)     
	 v_esc(39)       Phi0(40)        n_gas(41)       M_gas(42)       lambda_gas(43)  Lx_gas(44)      Ly_gas(45)      Lz_gas(46)      
	a_gas(47)       Eax_gas(48)     Eay_gas(49)     Eaz_gas(50)     b_gas(51)       Ebx_gas(52)     Eby_gas(53)     Ebz_gas(54)     
	 c_gas(55)       Ecx_gas(56)     Ecy_gas(57)     Ecz_gas(58)     Ekin_gas(59)    Epot_gas(60)    lambdaE_gas(61) 
	 n_star(62)      M_star(63)      lambda_star(64) Lx_star(65)     Ly_star(66)     Lz_star(67)     a_star(68)      
	Eax_star(69)    Eay_star(70)    Eaz_star(71)    b_star(72)      Ebx_star(73)    Eby_star(74)    Ebz_star(75)    c_star(76)      
	 Ecx_star(77)    Ecy_star(78)    Ecz_star(79)    Ekin_star(80)   Epot_star(81)   lambdaE_star(82) */

	/* Col:		    1   2    3  4  5  6  7  8  9 10 11 */
	sscanf(lineRead, "%d  %d  %d %f %d \
			  %f   %f   %f %f %f %f \
			  %f   %f   %f %f %f %f %f %f %f %f \
			  %f   %f   %f \
			  %f   %f   %f %f %f %f %f %f %f %f \
			  %f   %f   %f %f %f %f %f %f %f %f ",
 
			&dummyID, &halo->hostID, &halo->nSub, &halo->mTot, &tmpNpart, 
			&halo->X[0], &halo->X[1], &halo->X[2], &halo->V[0], &halo->V[1], &halo->V[2], 				// 11
			&halo->rVir, &dummy, &halo->rsNFW, &dummy, &dummy, &halo->vMax, &dummy, &halo->sigV, &halo->lambda, &dummy, // 21
			&halo->L[0], &halo->L[1], &halo->L[2],									// 24
			&dummy, &dummy, &dummy, &dummy,   &dummy, &dummy, &dummy, &dummy, &dummy, &dummy,			// 34
			&dummy, &dummy, &dummy, &halo->fMhires, &dummy, &dummy, &dummy, &dummy, &halo->cNFW, &dummy);		 // 44

	/* Particle numbers were not allocated correctly sometimes, so let's reset them carefully */
	nGas = 0; nStar = 0;

#ifdef IDADD
	if (dummyID < idADD)
		dummyID += dummyID;
#endif

	halo->ID = dummyID;
	halo->nPart[0] = nGas;
	halo->nPart[1] = tmpNpart - nGas - nStar;

	if (nPTypes > 2)
	{
		halo->nPart[2] = 0;
		
		if (nPTypes > 3)
		{
			halo->nPart[3] = 0;
		
			if (nPTypes > 4)
			{
				halo->nPart[4] = nStar;
					
					if (nPTypes > 5)
					{
						halo->nPart[5] = 0;
					}
			}
		}
	}

	/* Compute max velocity and sub box edges while reading the halo file ---> this is used to compute the buffer zones */
	vHalo = VectorModule(halo->V);

	if (vHalo > locVmax)
		locVmax = vHalo;
};
#endif



void IOSettings::WriteTree(int iThisCat)
{
	string outName;
        string strCpu = to_string(locTask);
	int orphan, iType = 0;

	iType = 1;

        for (int iC = iThisCat-1; iC < iThisCat; iC++)
        {
		const char *strFnm;	
		strFnm = strSnaps[iC].c_str();
		outName = pathOutput + outPrefix + strFnm + "." + strCpu + "." + outSuffix;

		ofstream fileOut;
		fileOut.open(outName);

                if (locTask == 0)
                        cout << "Printing trees to file " << outName << endl;

		if (locTask == 0)
		{
			fileOut << "# ID host(1)   N particles host(2)   Num. progenitors(3)  Orphan[0=no, 1=yes](4)" << endl;
			fileOut << "# Total particles (1)   ID progenitor(2)   Particles in common (3)" << endl;
		} 

                for (int iM = 0; iM < locCleanTrees[iC].size(); iM++)
                {
			MergerTree thisTree = locCleanTrees[iC][iM];

			if (thisTree.isOrphan)
				orphan = 1;	
			else
				orphan = 0;

			int nTotPt = 0;
			for (int iA = 0; iA < nPTypes; iA++)
				nTotPt += thisTree.mainHalo.nPart[iA];

			fileOut << thisTree.mainHalo.ID 	<< " " 
				<< nTotPt 			<< " " 
				<< thisTree.idProgenitor.size() << " "
				<< orphan << endl;

                        for (int iP = 0; iP < thisTree.idProgenitor.size(); iP++)
			{
				Halo progHalo = thisTree.progHalo[iP];

				int nTotPt = 0, nTotComm = 0;
		
				for (int iA = 0; iA < nPTypes; iA++)
					nTotPt += progHalo.nPart[iA];

				for (int iA = 0; iA < nPTypes; iA++)
					nTotComm += thisTree.nCommon[iA][iP];

				fileOut	<< nTotPt 		<< " " 
                                	<< progHalo.ID		<< " "
					<< nTotComm	 	<< endl;
			}
                }	// loop on merger tree
		
		fileOut.close();
        }	// loop on iCatalog
};


void IOSettings::WriteLog(int iNum, float time)
{
	/* Number of times stored in the log file */
	int nLogStep = 6;

	if (iNum == 0)
	{
        	string strCpu = to_string(totTask);
		string strChu = to_string(nChunks);

#ifdef ZOOM
		strChu += ".zoom";
#endif

		outLogName = pathOutput + "timing_log." + strCpu + "." + strChu + ".txt";
		fileLogOut.open(outLogName);
		fileLogOut << "# ReadFile (1) Communication (2) ForwardTree (3) BackwardTree (4) SyncBuffer(5) Memory (6)" << endl;
	} else if (iNum > 0) {
		logTime.push_back(time);
	
		if (logTime.size() == nLogStep) 
		{
			for (auto const& thisTime : logTime)
				fileLogOut << setw(9) << thisTime << "\t";
			
			/* End and clean the line */
			fileLogOut << endl;
			logTime.clear();
		}
	} else if (iNum == -1) {
		fileLogOut.close();
	}
};



tk::spline IOSettings::ReadA()
{
	double a = 0.0, t = 0.0, GYr = 0.005;
	vector <double> as, ts;
	int aStep = 0;
	string lineIn;
	
	tk::spline a2t; 

	ifstream fileIn(pathA);

	if (!fileIn.good())
	{
		cout << "File: " << pathA << " not found on task=" << locTask << endl;
	} else {
		if (locTask == 0)
       			cout << "Reading tree file: " << pathA << endl;
	}	

	while (getline(fileIn, lineIn))
	{
		const char *lineRead = lineIn.c_str();
		sscanf(lineRead, "%lf", &a);
		aStep++;

		t = aStep * GYr;
	
		ts.push_back(t);
		as.push_back(a);
	}

	a2t.set_points(as, ts);

	//cout << locTask << ") a2t test = " << a2t(0.9999) << endl;
	//cout << locTask << ") a2t test = " << a2t(0.5) << endl;

	if (locTask == 0)
		cout << "Read file: " << pathA << " with " << as.size() << " lines." << endl;

	return a2t;
};


tk::spline IOSettings::ReadPk()
{
	double k = 0.0, pk = 0.0;
	vector <double> ks, pks;
	string lineIn;
	
	tk::spline kPk; 

	ifstream fileIn(pathPk);

	if (!fileIn.good())
	{
		cout << "File: " << pathPk << " not found on task=" << locTask << endl;
	} else {
		if (locTask == 0)
       			cout << "Reading tree file: " << pathPk << endl;
	}	

	while (getline(fileIn, lineIn))
	{
		const char *lineRead = lineIn.c_str();
		sscanf(lineRead, "%lf  %lf", &k, &pk);

		ks.push_back(k);
		pks.push_back(pk);
	}

	kPk.set_points(ks, pks);

	if (locTask == 0)
		cout << "Read file: " << pathPk << " with " << ks.size() << " lines." << endl;


	//cout << locTask << ") Pk test = " << kPk(0.1) << endl;

	return kPk;

};


void IOSettings::WriteSmoothTrees()
{};
