#!/usr/bin/python3

from mtreelib.sqllib import *
from mtreelib.mtree import *
from mtreelib.read_mtree import *
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

''' 
	test_mtree.py
	This is a basic script to test the python post processing routines
'''

print('Testing Merger Tree python post-processing scripts')

nSnaps = 54	
nSteps = 4
nChunk = 1

thisCode = '00_06'


#baseTreeMCPP = '/home/edoardo/devel/MetroC++/output/fullbox_01_'
#baseTreeMCPP = '/z/carlesi/CLUES/MetroC++/output/HESTIA/2048/GAL_FOR/00_06/hestia_2048_00_06_'
baseTreeMCPP = '/home/eduardo/CLUES/MetroC++/output/2048/'+thisCode+'/00/out_'
suffTreeMCPP = 'mtree'

thisDb = 'trees_' + thisCode + '.db'
newSql = SQL_IO(thisDb, nSteps)

try:
	newSql.halo_table()
except ValueError:
	print("Halo table already exists in SQL database.")

readFiles = ReadSettings(baseTreeMCPP, suffTreeMCPP, nChunk, nSnaps, nSteps)
allTrees = readFiles.read_trees()

for thisTree in allTrees:
	[tmp_m, tmp_id] = thisTree.get_mass_id()
	newSql.insert_tree(tmp_id[0], tmp_m)

newSql.close()
print('Done.')

