#!/bin/bash

#mpiexec='mpiexec'		# Laptop test
#mpiexec='/z/eduardo/CLUES/libs/bin/mpiexec'
mpiexec='/home/eduardo/CLUES/libs/bin/mpiexec'

cd ..; make; cd -
rm ../tmp/*.tmp

cd .. ; cd bin/ ; mv MetroCPP MetroCPP-zoom ; cd ../scripts 

#$mpiexec -n 12 ../bin/MetroCPP ../config/fullbox_leibniz.cfg	
#$mpiexec -n 4 ../bin/MetroCPP ../config/lgf_tmp.cfg
$mpiexec -n 1 ../bin/MetroCPP-zoom ../config/lgf_2048.cfg
#$mpiexec -n 4 ../bin/MetroCPP ../config/fullbox_laptop.cfg	# Laptop setting
