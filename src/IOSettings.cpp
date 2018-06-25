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
				iLine++;
			}

		} else if (iLine == 1) {

	                sscanf(lineRead, "%u %llu", &nPartHalo, &locHaloID);
			iLine++;

		} else {

			locParts[iLocHalos][iLocParts].ReadLineAHF(lineRead);
			totLocParts++;
			iLocParts++;

			if (iLocParts == locHalos[iLocHalos].nPart)
			{
				iLocParts = 0;
				iLocHalos++;
				iLine = 1;	// Reset some indicators
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

	unsigned int iLocHalos = 0, nPartHalo = 0;
	
	nLocHalos = NumLines(urlHalo);	
	
	locParts.resize(nLocHalos); 
	locHalos.resize(nLocHalos); 

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
			locHalos[iLocHalos].ReadLineAHF(lineRead);
			nPartHalo = locHalos[iLocHalos].nPart;
			locParts[iLocHalos].resize(nPartHalo);
			locPartsSize += nPartHalo * sizePart;
			nLocParts += nPartHalo; 
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

