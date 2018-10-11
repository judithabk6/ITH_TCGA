#!/bin/bash

i=1
for PARAM in cancer_loc pythonenv
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done


source $pythonenv/bin/activate

./src/get_SNV_filters_and_data_exploration.py $cancer_loc
