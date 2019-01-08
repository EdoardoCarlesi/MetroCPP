#!/usr/bin/python3

from mtreelib.mtree import *
from mtreelib.read_mtree import *
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

''' 
	test_mtree.py
	This is a basic script to test the python post processing routines
'''

print('Testing Merger Tree python post-processing scripts')
	
nSteps = 4
nChunk = 4

baseTreeMCPP = '/home/edoardo/devel/MetroC++/output/fullbox_01_'
suffTreeMCPP = 'mtree'

readFiles = ReadSettings(baseTreeMCPP, suffTreeMCPP, nChunk, nSteps)
allTrees = readFiles.read_trees()


#for ih in range(0, 10):
#	allTrees[ih].print_mass_id()


print('Done.')

