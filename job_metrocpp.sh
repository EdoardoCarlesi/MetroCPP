#!/bin/bash

#PBS -q debugn
#PBS -r n
#PBS -l nodes=2:ppn=12
#PBS -l walltime=01:45:00
#PBS -u eduardo
#PBS -m abe
#PBS -N metro_cpp

#Your work dir is that folder where you submit the job
WDIR=/home/eduardo/CLUES/MetroCPP/

echo cd $WDIR
cd $WDIR

mcpp=bin/MetroCPP
cfgf=config/fullbox_1024.cfg

module load compilers/gcc/6.1.0

echo mpirun ./$mcpp $cfgf
mpirun ./$mcpp $cfgf
