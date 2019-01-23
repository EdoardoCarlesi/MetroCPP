#!/opt/rh/rh-python36/root/usr/bin/python3
#########################/usr/bin/python3

from mtreelib.sqllib import *
from mtreelib.mtree import *
from mtreelib.read_mtree import *
import timeit
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

''' 
	loop_treedir.py
	Extract trees and store them into a database
'''

print('Testing Merger Tree python post-processing scripts')

baseTreeMCPP = '/z/carlesi/STORE/LGF/trees/metrocpp/'
suffTreeMCPP = 'mtree'

nSnaps = 54	
nSteps = 54
nChunk = 1

iSeedIni = 0
iSeedEnd = 1
gSeedIni = 10
gSeedEnd = 20

# Save all the extracted trees into a database
thisDb = baseTreeMCPP + 'LGF_all_trees.db'

# Initialize the database and begin transaction. This avoids to commit at every step and speeds up the program
newSql = SQL_IO(thisDb, nSteps)
newSql.halo_table()
newSql.cursor.execute('BEGIN TRANSACTION')

# Loop on the large scale seed codes
for iSeed in range(iSeedIni, iSeedEnd):
	iSeedStr = '%02d' % iSeed

	# Loop on the small scale seed codes
	for gSeed in range(gSeedIni, gSeedEnd):
		gSeedStr = '%02d' % gSeed
		
		thisSubDir = iSeedStr + '_' + gSeedStr
		thisTreePath = baseTreeMCPP + thisSubDir
		rootFile = thisSubDir + '/lgf_' + thisSubDir + '_'	
		testFile = thisSubDir + '/lgf_' + thisSubDir + '_' + '%03d' % nSteps + suffTreeMCPP
	
		# If this path exists then extract the merger tree therein
		if os.path.isdir(thisTreePath) and os.path.isfile(thisTreePath+testFile):
			readFiles = ReadSettings(baseTreeMCPP+rootFile, suffTreeMCPP, nChunk, nSnaps, nSteps)
			allTrees = readFiles.read_trees()
		
			start = timeit.default_timer()
			print('Adding %s to database %s.' % (thisSubDir, thisDb))

			# Store all the trees inside a database
			for thisTree in allTrees:
				[tmp_m, tmp_id] = thisTree.get_mass_id()
				newSql.insert_tree(tmp_id[0], thisSubDir, tmp_m, tmp_id)

			end = timeit.default_timer()
			print('Inserted %d trees in %f seconds.' % (len(allTrees), end - start))

# Now commit to the database and close
newSql.cursor.execute('COMMIT')
newSql.close()
