#ifndef IOSETTINGS_H
#define IOSETTINGS_H

#include <vector>
#include <string>

using namespace std;

class IOSettings {

public:
	IOSettings();
	~IOSettings();

	int nLinesHeader = 1;

	int iniStep;
	int endStep;
	int totStep;
		
	int nFilesOnTask;

	// Catalog numbers, redshifts and expansion factor values
	vector<int> numSnaps;
	vector<string> strSnaps;
	vector<float> redShift;
	vector<float> aFactors;

	// These strings contain all the paths to all the halo files to be read
	// They are vector of vectors as in principle a task may read from different files
	vector <vector<string>> haloFiles;
	vector <vector<string>> partFiles;

	string haloSuffix;
	string partSuffix;
	string haloPrefix;
	string partPrefix;

	int nCat;		// Number of catalogs
	string catFormat;
	string thisPath;	// Where the program is installed
	string pathInput;	// Folder containing all the halo/particle catalogs	//TODO allow for different folders
	string pathOutput;	// Where to dump all the output data

	// Temporary storage for test URLs
	string urlTestFileHalo; 
	string urlTestFilePart; 
	
	void Init(void);
	void DistributeFilesAmongTasks(void);

	// These function rely on bash scripts to determine the halo catalogs properties
	void FindCatN();
	void FindCatZ();
	void FindCatID();

	void ReadConfigFile(string);

	void ReadParticles();
	void ReadHalos();
	void WriteTrees();

private:
	// TODO : AHF-halo catalog format only is assumed!
	// These scripts are used by the program to find:
	string findNsh = "scripts/find_n.sh";	// Number of catalogs
	string findZsh = "scripts/find_z.sh";	// Redshift of snapshots
	string findIDsh = "scripts/find_id.sh";	// ID number of catalogs
	string tmpIdOut = "tmp/output_id.tmp";
	string tmpZOut = "tmp/output_z.tmp";
	string tmpNOut = "tmp/output_n.tmp";
};




#endif
