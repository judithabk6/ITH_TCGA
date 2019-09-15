#!/bin/bash

i=1
for PARAM in url
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done

cd /data/users/jabecass/ITHeterogene/rerun_tcga_het_project/single_cell_dataset/CRC_leung/raw_data
wget $url