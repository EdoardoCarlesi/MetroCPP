#!/bin/bash

mpiexec='/home/eduardo/CLUES/libs/bin/mpiexec'

$mpiexec -n 4 ../bin/MetroCPP ../config/test_zoom.cfg
