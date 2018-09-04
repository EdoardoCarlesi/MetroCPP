#!/bin/bash
# Find the redshift of a given halo catalog

path=$1
format=$2
zoom=$3
cpus='0000'
cd $path

if [ "$zoom" == "true" ]
then
 	cpus='.'
else
	cpus='.0000.'
fi

n_snaps=`ls *$cpus*$format | wc -l`

echo $n_snaps
