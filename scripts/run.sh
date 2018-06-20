#!/bin/bash

#make clean; make
n_cpus=4

rm bin/*
make
mpirun -n $n_cpus ./bin/PTrees
