#ifndef IOSETTINGS_H
#define IOSETTINGS_H

#include <vector>
#include <string>
#include "Halo.h"

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
	void SetCosmology();

	/* Read input */
	void ReadPk();
	void ReadA();
	void ReadLineAHF(const char *, Halo *);
	void ReadParticles();
	void ReadHalos();
	void ReadTrees();

	/* Write output */
	void WriteTrees();
	void WriteSmoothTrees();

private:
	void InitFromCfgFile(vector<string>);

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
