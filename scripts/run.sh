#!/bin/bash

#make clean; make
n_cpus=12

rm bin/*
make
mpirun -n $n_cpus ../bin/MetroCPP ../config/fullbox_1024.cfg
#mpirun -n $n_cpus ./bin/MetroCPP config/config_template.cfg
