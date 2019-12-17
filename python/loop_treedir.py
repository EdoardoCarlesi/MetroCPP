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

#resolution='2048'
resolution='4096'
# This is assuming an input .mtree file structure like: hestia_RESOLUTION_XX_XX.0.mtree
#baseTreeMCPP = '/z/carlesi/STORE/LGF/trees/metrocpp/'
#baseTreeMCPP = '/home/eduardo/CLUES/DATA/HESTIA/8192/trees/'
#baseTreeMCPP = '/z/carlesi/CLUES/MetroC++/output/HESTIA/4096/'; baseRootFile = 'hestia_4096_'
baseTreeMCPP = '/home/eduardo/CLUES/DATA/trees/' + resolution + '/'; baseRootFile = 'out_'
suffTreeMCPP = 'mtree'

<<<<<<< HEAD
nSnaps = 127
nSteps = 126
=======
nSteps = 54

if resolution == '2048':
    nSnaps = 54
if resolution == '4096':
    nSnaps = 55

>>>>>>> ee8a49a50d9b46059c392e41d6993ee37b515e3f
nChunk = 1

# Print to file all halos with a z=0 particle content above this value
partThreshold = 20

iSeedIni = 0; iSeedEnd = 60; gSeedIni = 0; gSeedEnd = 20; 
#iSeedIni = 17; iSeedEnd = 38; gSeedIni = 11; gSeedEnd = 12; 
#iSeedIni = 18; iSeedEnd = 80; gSeedIni = 0; gSeedEnd = 20;
#iSeedIni = 37; iSeedEnd = 38; gSeedIni = 0; gSeedEnd = 10;

nSeedIni = 0; nSeedEnd = 10

# All the trees will be stored inside this database
thisDb = baseTreeMCPP + 'hestia_trees_' + resolution + '.db'

# Initialize the database and begin transaction. This avoids to commit at every step and speeds up the program
newSql = SQL_IO(thisDb, nSteps)
newSql.halo_table()
newSql.cursor.execute('BEGIN TRANSACTION')

# Loop on the large scale seed codes
for iSeed in range(iSeedIni, iSeedEnd):
    iSeedStr = '%02d' % iSeed 
#    iSeedStr = iSeedStr + '_11'

    # Loop on the small scale seed codes
    for gSeed in range(gSeedIni, gSeedEnd):
        gSeedStr = '%02d' % gSeed

        # Loop on the small scale seed codes
        for nSeed in range(nSeedIni, nSeedEnd):
            nSeedStr = '%02d' % nSeed

            thisSubDir = iSeedStr + '_' + gSeedStr + '/' + nSeedStr
            thisTreePath = baseTreeMCPP + thisSubDir + '/'

            # General file structure
            #rootFile = thisTreePath + baseRootFile + thisSubDir + '_'
            rootFile = thisTreePath + baseRootFile

            # We need to set a default test file to check whether the directory contains some actual stuff
            testFile = rootFile + '%03d' % nSnaps + '.0.' + suffTreeMCPP

            #print(rootFile, testFile)

            #if os.path.isdir(thisTreePath):
            #    print(thisTreePath)

            # If this path exists then extract the merger tree therein
            if os.path.isdir(thisTreePath) and os.path.isfile(testFile):
                readFiles = ReadSettings(rootFile, suffTreeMCPP, nChunk, nSnaps, nSteps)
                allTrees = readFiles.read_trees()

                start = timeit.default_timer()
                print('Adding %s to database %s.' % (thisSubDir, thisDb))

                # Store all the trees inside a database
                for thisTree in allTrees:
                    [tmp_m, tmp_id] = thisTree.get_mass_id()

                    # Save all the trees inside the DB
                    newSql.insert_tree(tmp_id[0], thisSubDir, tmp_m, tmp_id)

                    # Print to ASCII File
                    if tmp_m[0] > partThreshold:
                        this_file = thisTreePath + 'halo_' + str(tmp_id[0]) + '.ids'
                        thisTree.dump_to_file_id(this_file)

                end = timeit.default_timer()
                print('Inserted %d trees in %f seconds.' % (len(allTrees), end - start))

            else:
                'Do nothing'
            #print('Files not found: ', rootFile, testFile)

# Now commit to the database and close
newSql.cursor.execute('COMMIT')
newSql.close()
