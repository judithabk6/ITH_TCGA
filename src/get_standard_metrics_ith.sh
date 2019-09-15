#!/bin/bash

i=1
for PARAM in patient
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done


source /data/users/jabecass/dl_tools_centos/miniconda3_activate
source activate my_python36_centos

./src/get_standard_metrics_ith.py $patient
