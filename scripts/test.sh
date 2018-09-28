#!/bin/bash

mpiexec='/home/eduardo/CLUES/libs/bin/mpiexec'

#$mpiexec -n 4 ../bin/MetroCPP ../config/zoom_2048.cfg
$mpiexec -n 2 ../bin/MetroCPP ../config/test_zoom.cfg
