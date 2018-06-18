#!/bin/bash

make clean; make
mpirun -n 2 ./bin/PTrees
