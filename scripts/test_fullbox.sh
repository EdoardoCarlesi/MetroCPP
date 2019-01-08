#!/bin/bash

mpiexec='mpiexec'		# Laptop test
#mpiexec='/z/eduardo/CLUES/libs/bin/mpiexec'

cd ..; make; cd -
rm ../tmp/*.tmp

<<<<<<< HEAD
#$mpiexec -n 4 ../bin/MetroCPP ../config/fullbox_1024.cfg	# Neruda
#$mpiexec -n 4 ../bin/MetroCPP ../config/fullbox_laptop.cfg	# Laptop setting
$mpiexec -n 4 ../bin/MetroCPP ../config/fullbox_leibniz.cfg	
#$mpiexec -n 4 ../bin/MetroCPP ../config/fullbox_1024.cfg
#$mpiexec -n 4 ../bin/MetroCPP ../config/fullbox_laptop.cfg	# Laptop setting
=======
#$mpiexec -n 4 ../bin/MetroCPP ../config/fullbox_1024.cfg
$mpiexec -n 4 ../bin/MetroCPP ../config/fullbox_laptop.cfg	# Laptop setting
>>>>>>> 8defbc797aa01d41fb6cbcf33b69eecf89c9e713
