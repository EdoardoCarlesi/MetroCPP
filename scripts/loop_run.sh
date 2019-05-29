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

#00_06  01_12  09_18  17_10  17_13  34_13  37_11  45_17  55_02  62_14  64_14
RunNum[4]="01_12"
RunNum[5]="09_18"
RunNum[6]="17_10"
RunNum[7]="17_13"
RunNum[8]="37_11"
RunNum[9]="62_14"
RunNum[10]="64_14"

resolution='2048'

ncpu=1

iniRun=4
totRun=11

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
		sed 's/RUNNUM/'${Run}'/g' < $TemplCfg | sed 's/RESOLUTION/'$resolution'/g' | sed 's/SUBRUN/'${SubRun}'/g' &> $CfgTmp

		mkdir $MetroOut'/'$Run'/'$SubRun
		echo 'Sub-simulation number:' $CfgTmp
		echo mpirun -n $ncpu $MetroExe $CfgTmp
		mpirun -n $ncpu $MetroExe $CfgTmp
		#head $CfgTmp
	done
done
