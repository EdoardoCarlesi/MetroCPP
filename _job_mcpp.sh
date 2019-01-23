#!/bin/bash

#PBS -q QUEUE
#PBS -r n
#PBS -l nodes=NODES:ppn=PPN
#PBS -l walltime=HOURS:MINUTES:00
#PBS -u eduardo
#PBS -m abe
#PBS -N mcppRUNCODE

#Your work dir is that folder where you submit the job
WDIR=/home/eduardo/CLUES/MetroCPP/

echo cd $WDIR
cd $WDIR

runCode=RUNCODE

mcpp=bin/MetroCPP
base_cfg=config/fullbox_template_leibniz.cfg
tmp_cfg=config/tmp_fullbox_RUNCODE.cfg

module load compilers/gcc/6.1.0

#source /opt/env/intel-2017.sh
#module load mpi/mvapich2-1.9
#module load compilers/gcc/4.8.3   
#module load compilers/gcc/8.1.0
#module load compilers/gcc/4.8.3
#module load mpi/openmpi-2.0.2-gcc-4.8.3

echo mpirun ./$mcpp $tmp_cfg
mpirun ./$mcpp $tmp_cfg &> mcppRUNCODE_${PBS_JOBID}.out
