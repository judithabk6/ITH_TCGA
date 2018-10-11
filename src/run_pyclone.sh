#!/bin/bash

i=1
for PARAM in patient folder pythonenv
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done


source $pythonenv/bin/activate

purity=`cat results/$patient/pyclone/ascat_purity.txt`
PyClone run_analysis_pipeline --in_files results/$patient/pyclone/$folder/input.tsv --working_dir results/$patient/pyclone/$folder --tumour_contents $purity 

