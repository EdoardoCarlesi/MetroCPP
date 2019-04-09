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

#baseTreeMCPP = '/z/carlesi/STORE/LGF/trees/metrocpp/'
#baseTreeMCPP = '/home/eduardo/CLUES/DATA/LGF/trees/'
#baseTreeMCPP = '/home/eduardo/CLUES/DATA/HESTIA/8192/trees/'
baseTreeMCPP = '/z/carlesi/CLUES/MetroC++/output/HESTIA/4096/'
#suffTreeMCPP = 'mtree'
suffTreeMCPP = 'mtree'

nSnaps = 127
nSteps = 127
nChunk = 1

#iSeedIni = 17; iSeedEnd = 18; gSeedIni = 11; gSeedEnd = 12; thisRun = '17_11'
#iSeedIni = 37; iSeedEnd = 38; gSeedIni = 11; gSeedEnd = 12; thisRun = '37_11'
thisDb = baseTreeMCPP + 'hestia_trees_' + thisRun + '.db'


iSeedIni = 0; iSeedEnd = 2; gSeedIni = 0; gSeedEnd = 20; thisRun = '17_11'
thisDb = baseTreeMCPP + 'hestia_trees_4096.db'

# Save all the extracted trees into a database
#thisDb = baseTreeMCPP + 'lgf_n500_trees.db'
#thisDb = baseTreeMCPP + 'hestia_trees_test.db'

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
		#rootFile = thisTreePath + '/lgf_' + thisSubDir + '_'	
		rootFile = thisTreePath + '/hestia_8192_' + thisSubDir + '_'	
		testFile = rootFile + '%03d' % nSnaps + '.0.' + suffTreeMCPP

		'''
		#thisTreePath = baseTreeMCPP + 'trees/00_10mpi'
		#rootFile = thisTreePath + '/lgf_00_10_'	
		#testFile = thisTreePath + '/lgf_00_10_' + '%03d' % nSnaps + '.0.' + suffTreeMCPP
	
		thisTreePath = baseTreeMCPP + thisSubDir + '/'
		rootFile = thisTreePath	+ 'hestia_8192_17_11_'
		testFile = rootFile + '%03d' % nSnaps + '.0.' + suffTreeMCPP
		'''

		print(rootFile, testFile)

		# If this path exists then extract the merger tree therein
		if os.path.isdir(thisTreePath) and os.path.isfile(testFile):
			readFiles = ReadSettings(rootFile, suffTreeMCPP, nChunk, nSnaps, nSteps)
			allTrees = readFiles.read_trees()
		
			start = timeit.default_timer()
			print('Adding %s to database %s.' % (thisSubDir, thisDb))

			# Store all the trees inside a database
			for thisTree in allTrees:
                                [tmp_m, tmp_id] = thisTree.get_mass_id()

                                # Only save trees above a 500 particle threshold at z=0
                                #if tmp_m[0] > 500:
                                newSql.insert_tree(tmp_id[0], thisSubDir, tmp_m, tmp_id)
                                #print()
				    
                                # Print to ASCII File
				#thisTree.dump_to_file_mass_id()
				#print(tmp_id[0], tmp_m)

				#[tmp_m, tmp_id] = thisTree.get_mass_id()
                        #for thisTree in allTrees[0:10]:
                         #   print(thisTree.get_mass_id())

			end = timeit.default_timer()
			print('Inserted %d trees in %f seconds.' % (len(allTrees), end - start))

		else:
                        'DoNothing'
#			print(thisTreePath, testFile)

# Now commit to the database and close
newSql.cursor.execute('COMMIT')
newSql.close()
