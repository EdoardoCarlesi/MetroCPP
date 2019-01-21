#!/opt/rh/rh-python36/root/usr/bin/python3

from mtreelib.sqlio import *
from mtreelib.mtree import *
from mtreelib.read_mtree import *
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

''' 
	test_mtree.py
	This is a basic script to test the python post processing routines
'''

print('Testing Merger Tree python post-processing scripts')

nSnaps = 127	
nSteps = 3
nChunk = 1

#baseTreeMCPP = '/home/edoardo/devel/MetroC++/output/fullbox_01_'
baseTreeMCPP = '/z/carlesi/CLUES/MetroC++/output/HESTIA/2048/GAL_FOR/00_06/hestia_2048_00_06_'
suffTreeMCPP = 'mtree'

newSql = SQL_IO('test.db', nSteps)
#newSql.halo_table()

dummyPT1 = [10000, 9001, 8001]; dummyID1 = 12347937339393
dummyPT2 = [11000, 9000, 8000]; dummyID2 = 23347937339393
dummyPT3 = [21000, 1900, 1000]; dummyID3 = 13347937339393
dummyPT4 = [41000, 3900, 2800]; dummyID4 = 92347937339393

newSql.insert_tree(dummyID1, dummyPT1)
newSql.insert_tree(dummyID2, dummyPT2)
newSql.insert_tree(dummyID3, dummyPT3)
newSql.select_tree(dummyID1)

newSql.close()

'''
readFiles = ReadSettings(baseTreeMCPP, suffTreeMCPP, nChunk, nSnaps, nSteps)
allTrees = readFiles.read_trees()


for ih in range(0, 10):
	allTrees[ih].print_mass_id()
'''

print('Done.')

