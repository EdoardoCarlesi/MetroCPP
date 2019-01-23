#!/bin/bash

runCode='03'

base_cfg=config/fullbox_template_leibniz.cfg
tmp_cfg=config/tmp_fullbox_$runCode'.cfg'

base_job=job_mcpp.sh
tmp_job=tmp_mcpp.sh

queue='long40'
hours='04' minutes='30' nodes='4' ppn='10'

mkdir output/$runCode

sed 's/RUNNUM/'$runCode'/g' < $base_cfg &> $tmp_cfg

#sed 's/RUNCODE/'$runCode'/g' < $base_job | sed 's/HOURS/'$hours'/g' | sed 's/MINUTES/'$minutes'/g' | sed 's/NODES/'$nodes'/g' | \
#sed 's/PPN/'$ppn'/g' | sed 's/QUEUE/'$queue'/g' &> $tmp_job
#chmod +x $tmp_job
#qsub < $tmp_job
