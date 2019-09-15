#!/bin/bash

i=1
for PARAM in patient folder pythonenv
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done


source /data/users/jabecass/dl_tools_centos/miniconda3_activate
source activate my_python27_centos

if [[ $folder == *"absolute"* ]]; then
purity=`cat results/$patient/pyclone/absolute_purity.txt`
else
purity=`cat results/$patient/pyclone/ascat_purity.txt`
fi
PyClone run_analysis_pipeline --in_files results/$patient/pyclone/$folder/input.tsv --working_dir results/$patient/pyclone/$folder --tumour_contents $purity 

