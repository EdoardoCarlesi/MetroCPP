#!/bin/bash

this_run='01_12'
out_path='output/'

# Base path where the trees are located, file format and file base name
#this_dir='/store/clues/HESTIA/RE_SIMS/4096/GAL_FOR/'$this_run'/AHF_output/'	format='AHF'	file='HESTIA_100Mpc_4096_'$this_run'.127_halo_127'
this_dir='/z/carlesi/CLUES/MetroC++/output/4096/'$this_run'/'	format='MCPP'	file='halo_'$this_run'_127'


# Read some dwarf galaxy list and do the check ONLY for those

# Loop on all files and check for inconsistencies
for this_file in `ls $this_dir$file*`
do

# TEST FILES
#this_file='/store/clues/HESTIA/RE_SIMS/4096/GAL_FOR/'$this_run'/AHF_output/HESTIA_100Mpc_4096_'$this_run'.127_halo_127000000000008.dat'
#this_file='/z/carlesi/CLUES/MetroC++/output/4096/01_12/halo_01_12_127000000002350.full_tree'

n_lines=`wc -l $this_file`
echo $this_file $n_lines

# Now loop on all the trees 
python3 analyze_tree.py $this_file $format $out_path $this_run $n_lines

done
