#!/bin/bash
# Find the redshift of a given halo catalog. 
#			     snap cpu redshift
# Expecting format:  anyname_XXX.XXXX.zX.XXX.AHF_halos

path=$1
format=$2
zoom=$3

if [ "$zoom" == "true" ]
then
 	cpus='.'
else
	cpus='.0000.'
fi

tmp='ahf.tmp'
cd $path

if [ -f "$tmp_out" ]
then
	rm $tmp_out
fi

for f_ahf in `ls *$cpus*$format`
do
	echo $f_ahf > $tmp
	this_z=`grep -o z[0-9].[0-9][0-9][0-9] $tmp | sed 's/z//g ' `

	if [ "$this_z" == "" ]
	then
		this_z=`grep -o z[0-9][0-9].[0-9][0-9][0-9] $tmp | sed 's/z//g ' `
	fi

	echo $this_z
	#echo $this_z >> $tmp_out
done

if [ -f $tmp ]
then
	rm $tmp
fi

