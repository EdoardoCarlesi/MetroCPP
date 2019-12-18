from mtreelib.mtree import *
from mtreelib.read_mtree import *
import timeit 
import os 
import warnings 

warnings.simplefilter (action = 'ignore', category = FutureWarning) 

'''
        trees2halo.py
	Read the merger tree files and generate ASCII files for each halo above particle threshold
'''

print ('Testing Merger Tree python post-processing scripts')

#This is assuming an input .mtree file structure like: hestia_RESOLUTION_XX_XX.0.mtree
#baseTreeMCPP = '/z/carlesi/STORE/LGF/trees/metrocpp/'
suffTreeMCPP = 'mtree'

#There is a step less than the number of snapshots
nSnaps = 54
nSteps = 11

#In large box simulations with MPIO
nChunk = 1

#Print to file all halos with a z=0 particle content above this value
partThreshold = 2000

#thisTreePath = '/path/where/the/merger/trees/are/stored'
thisTreePath = '/z/carlesi/CLUES/MetroC++/output/' 
rootFile = thisTreePath + 'lgf_test_'

#If this path exists then extract the merger tree therein
#if os.path.isdir (thisTreePath) and os.path.isfile (testFile):
if os.path.isdir(thisTreePath):
	readFiles = ReadSettings(rootFile, suffTreeMCPP, nChunk, nSnaps, nSteps) 
	allTrees = readFiles.read_trees()
	start = timeit.default_timer()

#Store all the trees inside a database
for thisTree in allTrees:
	[tmp_m, tmp_id] = thisTree.get_mass_id()

	#Print to ASCII File
	if tmp_m[0] > partThreshold:
		this_file = thisTreePath + '/HALOS/halo_' + str (tmp_id[0]) + '.ids'
		thisTree.dump_to_file_id (this_file)

	end = timeit.default_timer()

else:
	print ('Path: ', thisTreePath, ' not found.')

print('MAHs of individual halos printed in %f seconds.' % (end - start))
