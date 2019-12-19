#!/bin/bash
# Find the redshift of a given halo catalog. 
#			     snap cpu redshift
# Expecting format:  prefix_XXX.XXXX.zX.XXX.AHF_halos

path=$1
format=$2
cpus=$3
split=$4 
prefix=$5

tmp=`pwd`'/ahf.tmp'
cd $path

if [ -f "$tmp_out" ]
then
	rm $tmp_out
fi

for f_ahf in `ls $prefix*$cpus*$format`
do
	echo $f_ahf > $tmp

if [ "$split" == "_" ]
then
	this_z=`grep -o _[0-9][0-9][0-9] $tmp | sed 's/_//g ' `
elif [ "$split" == ".z" ]
then
	this_z=`grep -o [0-9][0-9][0-9]'.z' $tmp | sed 's/\.z//g ' `
fi
	echo $this_z
done

if [ -f $tmp ]
then
	rm $tmp
fi

