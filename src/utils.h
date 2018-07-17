#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include <vector>
#include <algorithm>
#include <string>

#include "global_vars.h"
#include "utils.h"

/*
 * General scope functions and utilities
 */

unsigned int NumLines(const char *);

float VectorModule(float *);

void InitLocVariables(void);

void CleanMemory(int);

void ShiftHalosPartsGrids(void);

bool CompareHalos(int, int, int, int);

void FindProgenitors(int, int);

vector<string> SplitString(string, string);

vector<int> CommonParticles(vector<vector<unsigned long long int>>, 
vector<vector<unsigned long long int>>);


#endif
