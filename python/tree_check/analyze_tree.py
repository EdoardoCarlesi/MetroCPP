'''
	This python script analyzes a merger tree and looks for inconsistencies such as:
		- early mtree truncation
		- big discontinuity in the tree
'''

import numpy as np
from read_tree import *
import sys

# These data are relative to the input file
f_name = sys.argv[1]
f_type = sys.argv[2]

# Output error log
f_out_path = sys.argv[3]
f_out_code = sys.argv[4]

# Number of steps in the tree file
n_step = sys.argv[5]

f_err_log = f_out_path + 'treecheck_' + f_type + '_' + f_out_code + '_err.log'
f_err = open(f_err_log, 'a')

# This parameter defines a discontinous / suspicious MF
ratio_thr = 5.0

# If a halo with at least this many particles is truncated, send a warning sign!
part_thr = 1000.0

# Step threshold, signals when a halo above a certain number of particles is tracked for less than this number of steps
step_thr = 80

if f_type == 'AHF':
	this_tree = read_ahf_tree(f_name, n_step)

if f_type == 'MCPP':
	this_tree = read_mcpp_tree(f_name, n_step)

i_step = 0
parts_old = this_tree[0]

# Sanity check on the tree
for n_parts in this_tree:
	parts_new = n_parts

	p_ratio1 = float(parts_new - parts_old) / float(parts_old)
	p_ratio2 = float(parts_new - parts_old) / float(parts_new)
	p_r1_abs = np.sqrt(p_ratio1 * p_ratio1)
	p_r2_abs = np.sqrt(p_ratio2 * p_ratio2)

	# If the 
	i_step += 1
	if p_r1_abs > ratio_thr or p_r2_abs > ratio_thr: 
		f_err.write('Discontinuous MTree for file: %s\n ---> Line: %d, mRatio1: %f, mRatio2: %f\n' % (f_name, i_step, p_r1_abs, p_r2_abs))

	parts_old = parts_new

if n_parts > part_thr and i_step < step_thr:
	f_err.write('Truncation of MTree for file: %s\n ---> Line %d, %d particles.\n' % (f_name, i_step, int(parts_new)))

f_err.close()
