#ifndef IOSETTINGS_H
#define IOSETTINGS_H

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

	// Temporary storage for test URLs
	string urlTestFileHalo; 
	string urlTestFilePart; 
	
	string pathInput;
	string pathOutput;

	void Init(void);
	void DistributeFilesAmongTasks(void);

	void ReadParticles();
	void ReadHalos();
	void WriteTrees();
};




#endif
