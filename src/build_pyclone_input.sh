#!/bin/bash

i=1
for PARAM in patient folder filename
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done

source config.txt
source $python35_env_path/bin/activate

./src/build_pyclone_input.py $patient $folder $filename

