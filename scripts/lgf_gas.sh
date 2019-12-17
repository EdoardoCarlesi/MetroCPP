#!/bin/bash

MetroDir=$HOME/CLUES/MetroCPP/
TemplCfg=$MetroDir'config/lgf_gas.cfg'
CfgTmp=$MetroDir'config/lgf_tmp.cfg'
#MetroExe=$MetroDir'bin/MetroCPP'
MetroExe=$MetroDir'bin/MetroCPP-Zoom'
MetroOut=$MetroDir'output/37_11/'
MetroTmp=$MetroDir'tmp/'

nCpus=1

for (( iPath=0; iPath<1; iPath++ )) 
do

subDir=`printf %02d $iPath`

#mkdir $MetroOut'/'$subDir

		sed 's/SUBPATH/'${subDir}'/g' < $TemplCfg &> $CfgTmp

#		head -n 30 $CfgTmp
		echo 'Sub-simulation number:' $CfgTmp
		echo mpirun -n $ncpu $MetroExe $CfgTmp
		mpirun -n $nCpus $MetroExe $CfgTmp
done
