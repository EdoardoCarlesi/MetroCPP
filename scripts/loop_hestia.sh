#!/bin/bash

#resolution='4096'
resolution='2048'

MetroDir=/z/carlesi/CLUES/MetroC++/
TemplCfg=$MetroDir'config/hestia_template.cfg'
CfgTmp=$MetroDir'config/hestia_tmp.cfg'
MetroExe=$MetroDir'bin/MetroCPP'
MetroOut=$MetroDir'output/HESTIA/'$resolution'/'
MetroTmp=$MetroDir'tmp/'

HestiaBaseDir=/z/carlesi/HestiaNoam/RE_SIMS/

ThisDir=$HestiaBaseDir'/'$resolution'/GAL_FOR/'
#ThisDir=$HestiaBaseDir'/'$resolution'/DM_ONLY/'

iRun=0

for ThisRunNum in `ls $ThisDir | grep [0-9][0-9]_[0-9][0-9]`
do
#echo $ThisRunNum
RunNum[${iRun}]=$ThisRunNum iRun=$(expr "$iRun" + "1")
done

echo 'Total number of subfolders: ' $iRun

ncpu=1

iniRun=0
#totRun=1
totRun=$iRun

for (( iRun=$iniRun; iRun<$totRun; iRun++ ))
do
	Run=${RunNum[${iRun}]}
	echo 'Realisation: ' $Run
	
	
	mkdir $MetroOut'/'$Run

	echo rm $MetroTmp'/'*'.tmp'
	rm $MetroTmp'/'*'.tmp'
	sed 's/RESOLUTION/'$resolution'/g' < $TemplCfg | sed 's/HESTIARUN/'${Run}'/g' &> $CfgTmp

	mkdir $MetroOut'/'$Run'/'$SubRun
	echo 'Sub-simulation number:' $CfgTmp
	echo mpirun -n $ncpu $MetroExe $CfgTmp
	mpirun -n $ncpu $MetroExe $CfgTmp
done
