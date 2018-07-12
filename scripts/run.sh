#!/bin/bash

#make clean; make
n_cpus=4

rm bin/*
make
mpirun -n $n_cpus ./bin/MetroCPP config/config_template.cfg
