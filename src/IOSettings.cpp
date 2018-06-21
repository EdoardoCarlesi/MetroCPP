#include <string>
#include <sstream>
#include <fstream>

#include "IOSettings.h"
#include "general.h"

using namespace std;


IOSettings::IOSettings() 
{
};


IOSettings::~IOSettings() 
{
};


void IOSettings::Init()
{
	// Initialize all the main variables, checking for consistency
};


void IOSettings::DistributeFilesAmongTasks(void)
{
};


// Particle sizes have already been allocated in the ReadHalos() routines, do a safety check for the size
void IOSettings::ReadParticles(void)
{
	// TODO use read(buffer,size) to read quickly blocks of particles all at the same time 
	string strUrlPart = pathInput + urlTestFilePart;	// FIXME : this is only a test url for the moment
	const char *urlPart = strUrlPart.c_str();
	unsigned int iLocHalos = 0, iLocParts = 0, totLocParts = 0, iLine = 0, nPartHalo = 0, nFileHalos = 0;		
	unsigned long long int locHaloID;
	string lineIn;

	ifstream fileIn(urlPart);
		
	if (!fileIn.good())
		cout << "File: " << urlPart << " not found on task=" << locTask << endl;
	else 
	        cout << "Task=" << locTask << " is reading " << nLocParts << " particles from file: " << urlPart << endl;
	
	while (getline(fileIn, lineIn))
	{
		const char *lineRead = lineIn.c_str();		

		if (iLine == 0)
		{
	                sscanf(lineRead, "%d", &nFileHalos);
	
			if (nFileHalos != nLocHalos && iLine == 0)
			{
				cout << "WARNING on task=" << locTask << " expected: " << nLocHalos << " found: " <<  nFileHalos << endl;
				iLine++;
			} else {
				cout << "Particle file contains " << nLocHalos << " halos." << endl; 
				iLocHalos++;
				iLine++;
			}

		} else if (iLine == 1) {

	                sscanf(lineRead, "%u %llu", &nPartHalo, &locHaloID);
			//locParts[iLocHalos-1].resize(nPartHalo); // This way seems to be slightly slower
			iLine++;
			
		} else {

			Particle thisPart;
			thisPart.ReadLineAHF(lineRead);
			locParts[iLocHalos-1].push_back(thisPart); 
			//locParts[iLocHalos-1][iLocParts].ReadLineAHF(lineRead);
			totLocParts++;
			iLocParts++;

			if (iLocParts == locHalos[iLocHalos-1].nPart)
			{
				//locHalos[iLocHalos-1].Part = locParts[iLocHalos-1];
				iLine = 1;	// Reset some indicators
				iLocParts = 0;
			}
		} // else iLine not 0 or 1
	} // end while


};
 

/* Using AHF by default */
void IOSettings::ReadHalos()
{
	string strUrlHalo = pathInput + urlTestFileHalo;	// FIXME : this is only a test url for the moment
	const char *urlHalo = strUrlHalo.c_str();
	const char *lineHead = "#";
	string lineIn;

	unsigned int iLocHalos = 0;
	
	nLocHalos = NumLines(urlHalo);	
	
	locParts.resize(nLocHalos); 

	locPartsSize = 0;	// This will be increased while reading the file
	locHalosSize = sizeHalo * nLocHalos;

	ifstream fileIn(urlHalo);
	
	if (!fileIn.good())
		cout << "File: " << urlHalo << " not found on task=" << locTask << endl;
	else 
	        cout << "Task=" << locTask << " is reading " << nLocHalos << " halos from file: " << urlHalo << endl;

	while (getline(fileIn, lineIn))
	{
		const char *lineRead = lineIn.c_str();		
	
		if (lineRead[0] != lineHead[0])
		{
			Halo thisHalo; 
			thisHalo.ReadLineAHF(lineRead);
			locHalos.push_back(thisHalo);
			locPartsSize += locHalos[iLocHalos].nPart * sizePart;
			nLocParts += locHalos[iLocHalos].nPart;
			iLocHalos++;
		}
	}

	fileIn.close();

	cout << "On task=" << locTask << " " << locPartsSize/1024/1024 << " MB pt and " << locHalosSize/1024/1024 << " MB hl " << endl; 
	cout << "NHalos: " << locHalos.size() << " on task=" << locTask << endl;
};

   
void IOSettings::WriteTrees()
{

};

