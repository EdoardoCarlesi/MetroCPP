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

float* UnitVector(float *);

void InitLocVariables(void);

void CleanMemory(int);

void ShiftHalosPartsGrids(void);

vector<string> SplitString(string, string);

vector<int> SortIndexes(vector<float>);
#endif
