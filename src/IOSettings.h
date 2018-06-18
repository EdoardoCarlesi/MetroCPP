#ifndef IOSETTINGS_H
#define IOSETTINGS_H

#include <string>

using namespace std;

class IOSettings {

public:
	IOSettings();
	~IOSettings();

	int iniStep;
	int endStep;
	int totStep;

	string outputPath;
	string inputPath;

	void Init(void);
	void DistributeFilesAmongTasks(void);

	void ReadParticles();
	void ReadHalos();
	void WriteTrees();
};




#endif
