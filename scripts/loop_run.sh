#!/bin/bash

MetroDir=/home/eduardo/CLUES/MetroC++/
TemplCfg=$MetroDir'config/loop_template.cfg'
CfgTmp=$MetroDir'config/loop_tmp.cfg'
MetroExe=$MetroDir'bin/MetroCPP'
MetroOut=$MetroDir'output/2048/'
MetroTmp=$MetroDir'tmp/'

RunNum[0]="00_06"
RunNum[1]="45_17"
RunNum[2]="55_02"
RunNum[3]="34_13"

ncpu=1

iniRun=0
totRun=4

iniSub=0
totSub=10

for (( iRun=$iniRun; iRun<$totRun; iRun++ ))
do
	
	Run=${RunNum[${iRun}]}
	echo 'Realisation: ' $Run
	
	mkdir $MetroOut'/'$Run

	for (( iSub=$iniSub; iSub<$totSub; iSub++ ))
	do
		echo rm $MetroTmp'/'*'.tmp'
		rm $MetroTmp'/'*'.tmp'
		SubRun=`printf "%02d" ${iSub}`
		sed 's/RUNNUM/'${Run}'/g' < $TemplCfg | sed 's/SUBRUN/'${SubRun}'/g' &> $CfgTmp

		mkdir $MetroOut'/'$Run'/'$SubRun
		echo 'Sub-simulation number:' $CfgTmp
		echo mpirun -n $ncpu $MetroExe $CfgTmp
		mpirun -n $ncpu $MetroExe $CfgTmp
		#head $CfgTmp
	done
done
