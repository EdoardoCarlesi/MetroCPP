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

printAscii = True
#printAscii = False

#baseTreeMCPP = '/z/carlesi/STORE/LGF/trees/metrocpp/'
#baseTreeMCPP = '/home/eduardo/CLUES/DATA/HESTIA/8192/trees/'

# This is assuming an input .mtree file structure like: hestia_RESOLUTION_XX_XX.0.mtree
resolution = '4096'
#resolution = '8192'
baseTreeMCPP = '/z/carlesi/CLUES/MetroC++/output/HESTIA/'+resolution+'/'
baseRootFile = 'hestia_'+resolution+'_'
suffTreeMCPP = 'mtree'

nSnaps = 127
nSteps = 127
nChunk = 1

# Print to file all halos with a z=0 particle content above this value
partThreshold = 250

#iSeedIni = 0; iSeedEnd = 2; gSeedIni = 0; gSeedEnd = 20; 
iSeedIni = 17; iSeedEnd = 30; gSeedIni = 0; gSeedEnd = 20; 

# All the trees will be stored inside this database
thisDb = baseTreeMCPP + 'hestia_trees_' + resolution + '.db'

# The individual mtrees for each z=0 halo will be printed here
outPath = '/z/carlesi/CLUES/MetroC++/output/'+resolution+'/'

dwPath = '/store/erebos/nil/HESTIA/ANALYSIS/'+resolution+'_GAL_FOR/'

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
		thisTreePath = baseTreeMCPP + thisSubDir + '/'

		# General file structure
		rootFile = thisTreePath + baseRootFile + thisSubDir + '_'	
		
		# We need to set a default test file to check whether the directory contains some actual stuff
		testFile = rootFile + '%03d' % nSnaps + '.0.' + suffTreeMCPP

		#print(rootFile, testFile)

		# If this path exists then extract the merger tree therein
		if os.path.isdir(thisTreePath) and os.path.isfile(testFile):
			# List of dwarf galaxies
			dwList = dwPath + 'dwarfs_' + thisSubDir + '.txt'

			print(dwList)
			# Read the list of dwarf galaxies to be tracked
			dwIDs = read_dwarfs(dwList)	
	
			readFiles = ReadSettings(rootFile, suffTreeMCPP, nChunk, nSnaps, nSteps)
			allTrees = readFiles.read_trees()
		
			start = timeit.default_timer()
			print('Adding %s to database %s.' % (thisSubDir, thisDb))

			# Store all the trees inside a database
			for thisTree in allTrees:
				[tmp_m, tmp_id] = thisTree.get_mass_id()

				# TODO only save trees above a certain threshold!!!!
                                # Save all the trees inside the DB
				newSql.insert_tree(tmp_id[0], thisSubDir, tmp_m, tmp_id)
				idKey = str(tmp_id[0])

				# Print those halos whose IDs are in the Dwarf Galaxy list ONLY
				#if idKey in dwIDs.keys():
				if tmp_m[0] > partThreshold:
			#		print('Found key: ', idKey)
					thisOutPath = outPath + thisSubDir + '/'
					thisTree.dump_to_file_mass_id(thisOutPath, thisSubDir)
				#else:
				#	print('Could not find key: ', idKey)

			end = timeit.default_timer()
			print('Inserted %d trees in %f seconds.' % (len(allTrees), end - start))

		else:
			'Do nothing'
			#print('Files not found: ', rootFile, testFile)

# Now commit to the database and close
newSql.cursor.execute('COMMIT')
newSql.close()
