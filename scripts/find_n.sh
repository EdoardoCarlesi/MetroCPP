#!/bin/bash
# Find the redshift of a given halo catalog

path=$1
format=$2
cpus='0000'
cd $path

n_snaps=`ls *.$cpus*.$format | wc -l`

echo $n_snaps
