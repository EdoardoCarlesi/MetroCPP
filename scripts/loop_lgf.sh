#!/bin/bash

MetroDir=/z/carlesi/CLUES/MetroC++/
TemplCfg=$MetroDir'config/lgf_template.cfg'
CfgTmp=$MetroDir'config/lgf_tmp.cfg'
MetroExe=$MetroDir'bin/MetroCPP'
MetroOut=/z/carlesi/STORE/LGF/trees/metrocpp/
MetroTmp=$MetroDir'tmp/'
ThisDir=/z/carlesi/STORE/LGF/catalogs/AHF/

MinSize=1024

iRun=0

for ThisRunNum in `ls $ThisDir | grep [0-9][0-9]'_'[0-9][0-9]`
do
#echo $ThisRunNum
RunNum[${iRun}]=$ThisRunNum iRun=$(expr "$iRun" + "1")
done

echo 'Total number of subfolders: ' $iRun

rm -rf $MetroTmp'/'*'.tmp'	

ncpu=1
do_tree="true"
iniRun=0
totRun=1 #$iRun

for (( iRun=$iniRun; iRun<$totRun; iRun++ ))
do
	Run=${RunNum[${iRun}]}
	echo 'Realisation: ' $Run
	ThisOutDir=$MetroOut'/'$Run

		if [ -d $ThisOutDir ]; then
			thisSize=`du $ThisOutDir | awk '{print $1}'`
			
			#echo 'SIZE:' $thisSize

			if [ $thisSize -gt "$MinSize" ]; then
				do_tree='false'
				echo 'Skipping' $ThisOutDir ' of size: ' $thisSize ', do_tree: ' $do_tree
			else
				do_tree='true'
			fi			
		
		else			
			do_tree='true'
			mkdir $ThisOutDir
		fi

	is_tmp=0
	cd $MetroTmp
	
	case "$f" in *.tmp | *.tmp )
		is_tmp=`ls $MetroTmp/*.tmp | wc -l`
        	;;
	*)
	esac
	cd -

	if [ "$is_tmp" -gt "0" ]; then
		echo rm $MetroTmp'/'*'.tmp'
		rm $MetroTmp'/'*'.tmp'
	fi

	if [ "$do_tree" == "true" ]; then

		sed 's/HESTIARUN/'${Run}'/g' < $TemplCfg &> $CfgTmp
		echo 'Sub-simulation number:' $CfgTmp
		echo mpirun -n $ncpu $MetroExe $CfgTmp
		mpirun -n $ncpu $MetroExe $CfgTmp
	fi
done
