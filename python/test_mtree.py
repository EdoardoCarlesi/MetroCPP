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
nSteps = 54
nChunk = 1

thisCode = '00_06'

#baseTreeMCPP = '/home/edoardo/devel/MetroC++/output/fullbox_01_'
#baseTreeMCPP = '/z/carlesi/CLUES/MetroC++/output/HESTIA/2048/GAL_FOR/00_06/hestia_2048_00_06_'
baseTreeMCPP = '/home/eduardo/CLUES/MetroC++/output/2048/'+thisCode+'/00/out_'
suffTreeMCPP = 'mtree'

thisDb = 'trees_' + thisCode + '.db'
newSql = SQL_IO(thisDb, nSteps)

'''
# Fill the DB
newSql.halo_table()

readFiles = ReadSettings(baseTreeMCPP, suffTreeMCPP, nChunk, nSnaps, nSteps)
allTrees = readFiles.read_trees()

test_ids = []

#for thisTree in allTrees[0:10]:
for thisTree in allTrees:
	[tmp_m, tmp_id] = thisTree.get_mass_id()
	newSql.insert_tree(tmp_id[0], thisCode, tmp_m, tmp_id)
#	test_ids.append(tmp_id[0])


# Read from the DB
'''

#testID1='2941750477566971389'
#testID2='2938536608351584561'
testID3='5197802116663954028'

thisTree = newSql.get_full_mtree(testID3)
#thisTree.print_mass_id()
print(thisTree.norm_mass())


'''
for testID in test_ids:
	line = newSql.select_tree(testID)
	print(line[0])
	theseIDs = line[1].split()	

	for thisID in theseIDs:
		intID = (thisID.replace(",", ""))
		print(intID)
	#print(newSql.select_tree(testID))
	#print(line[1].split())
	#line.split()
'''

newSql.close()
print('Done.')
