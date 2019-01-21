#!/bin/bash
# Find the redshift of a given halo catalog. 
#			     snap cpu redshift
# Expecting format:  anyname_XXX.XXXX.zX.XXX.AHF_halos

path=$1
format=$2
zoom=$3
#cpus=''
cpus='0000'
tmp=`pwd`'/ahf.tmp'
cd $path

if [ "$zoom" == "true" ]
then
 	cpus='.'
else
	cpus='.0000.'
fi

if [ -f "$tmp_out" ]
then
	rm $tmp_out
fi

for f_ahf in `ls *$cpus*$format`
do
	echo $f_ahf > $tmp
	this_z=`grep -o _[0-9][0-9][0-9] $tmp | sed 's/_//g ' `
	#this_z=`grep -o [0-9][0-9][0-9]'.z' $tmp | sed 's/\.z//g ' `

	echo $this_z
done

if [ -f $tmp ]
then
	rm $tmp
fi

