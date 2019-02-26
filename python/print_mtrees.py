#!/usr/bin/python
####!/usr/share/Modules/tools/python3.6
####!/opt/rh/rh-python36/root/usr/bin/python3

import time
import pickle
from mtreelib.sqllib import *
from mtreelib.mtree import *
#import pandas as pd
#import subprocess as sbp
import os

box_size = 100000.0
base_path = '/z/carlesi/CLUES/MetroC++/output/HESTIA/8192/'
this_db='hestia_trees.db'

tot_files = 1
use_files = 1
n_steps = 127
n_lgs = 1

# Load the SQL database containing all the LGF trees
in_db = base_path + this_db
newSql = SQL_IO(in_db, n_steps)
columnReadID = 'allHaloIDs'
columnReadPT = 'allNumPart'

# Read the halo IDs at z=0
id_file_name='/store/erebos/nil/HESTIA/ANALYSIS/8192_GAL_FOR/dwarfs_17_11.txt'
id_file = open(id_file_name, 'r')
id_data = id_file.read().splitlines()

all_ids =[]
i_line = 0
iBroken = 0
iValid = 0

for line in id_data:
	line_id = line.split()
		
	if i_line > 0:
		all_ids.append(line_id[0])

	i_line +=1
#	print(line_id)

#for this_id in all_ids:
for this_id in all_ids[0:3]:
#	print(this_id, columnReadPT)
	this_tree = newSql.select_tree(this_id, columnReadPT)
#	this_ids = newSql.select_tree(this_id, columnReadID)
	valid_tree = np.where(this_tree > 0)

	this_file_name = 'halo_' + str(this_id) + '.full_mtree'
	this_file = open(this_file_name, 'w')


	if valid_tree[0].size > 40:
		iValid +=1 
		this_mtree = MergerTree(n_steps, this_tree)
	
		#print >> this_file, this_mtree.smooth_tree()
			#print(this_mtree.last_major_merger(True), this_mtree.last_major_merger(False))
	else:
		iBroken += 1

		
print('Found %d valid, %d broken trees.' % (iValid, iBroken))

newSql.close()
