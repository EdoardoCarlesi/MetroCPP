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


#ifndef IOSETTINGS_H
#define IOSETTINGS_H

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include "spline.h"
#include "Halo.h"
#include "Cosmology.h"

using namespace std;

class IOSettings {

public:
	IOSettings();
	~IOSettings();

	// Catalog numbers, redshifts and expansion factor values
	vector<string> strSnaps;
	vector<float> redShift;
	vector<float> aFactors;
	vector<int> numSnaps;

	// These strings contain all the paths to all the halo files to be read
	// They are vector of vectors as in principle a task may read from different files
	vector <vector<string>> haloFiles;
	vector <vector<string>> partFiles;

	string haloSuffix;
	string partSuffix;
	string haloPrefix;
	string partPrefix;
	string outPrefix;
	string outSuffix;

	string cpuString;	
	string splitString;

	string catFormat;
	string pathMetroCpp;	// Where the program is installed
	string pathInput;	// Folder containing all the halo/particle catalogs	//TODO allow for different folders
	string pathOutput;	// Where to dump all the output data
	string pathTree;

	string inputFormat;	// AHF, FoF, or anything else	// TODO only working with AHF for the moment

	void Init(void);
	void DistributeFilesAmongTasks(void);
	void CheckStatus(void);

	/* These function rely on bash scripts to determine the halo catalogs properties */
	void FindCatN();
	void FindCatZ();
	void FindCatID();

	/* Read and parse the configuration file */
	void ReadConfigFile(string);
	
	/* Cosmology functions */
	void SetCosmology(Cosmology*);

	/* Read input */
	void ReadLineAHF(const char *, Halo *);
	void ReadParticles();
	void ReadHalos();
	void ReadTrees();

	/* Write output */
	void WriteLog(int, float);
	void WriteTree(int);
	//void WriteTrees();
	void WriteSmoothTrees();

private:
	/* Log file properties */
	int iLogStep;
	string outLogName;
	ofstream fileLogOut;
	vector<float> logTime;

	void InitFromCfgFile(vector<string>);

	/* These functions read and interpolate from the right cosmological functions */
	tk::spline ReadPk();
	tk::spline ReadA();

	/* These variables will be set internally */
	string pathPk;
	string pathA;

	/* These scripts are used by the program to determine input file properties */
	string findNsh  = "/scripts/find_n.sh";		// Number of catalogs
	string findZsh  = "/scripts/find_z.sh";		// Redshift of snapshots
	string findIDsh = "/scripts/find_id.sh";	// ID number of catalogs

	/* Temporary output files */
	string tmpIdOut = "/tmp/output_id.tmp";
	string tmpZOut  = "/tmp/output_z.tmp";
	string tmpNOut  = "/tmp/output_n.tmp";

	/* These files will be read and used for interpolation in Cosmology */
	string dataPkPlanck = "/data/pk_planck.dat";
	string dataPkWMAP7  = "/data/pk_wmap7.dat";
	string dataAPlanck  = "/data/a2t_5Myr_planck.dat";
	string dataAWMAP7   = "/data/a2t_5Myr_wmap7.dat";
};




#endif
