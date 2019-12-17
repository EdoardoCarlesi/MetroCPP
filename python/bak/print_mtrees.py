import time
import pickle
from mtreelib.sqllib import *
from mtreelib.mtree import *
import os

resolution = '8192'
baseTreeMCPP = '/z/carlesi/CLUES/MetroC++/output/HESTIA/'+resolution+'/'
baseRootFile = 'hestia_'+resolution+'_'
suffTreeMCPP = 'mtree'
thisDb = baseTreeMCPP + 'hestia_trees_' + resolution + '.db'
nSteps = 1

thisSubDir = '17_11'

print('Loading DB: ', thisDb)

testID1 = 127000000000003
#         '127000000000359'
newSql = SQL_IO(thisDb, nSteps)
#newSql.halo_table()
#thisTree = newSql.get_full_mtree(testID1)
#thisTree = newSql.select_tree(testID1, 'allNumPart')
#thisTree = newSql.select_tree(testID1, 'allHaloIDs')
#thisTree.print_mass_id()
#print(thisTree) #.norm_mass())
print(newSql)

#thisTree.dump_to_file_mass_id(thisOutPath, thisSubDir)


