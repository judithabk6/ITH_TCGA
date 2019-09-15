#!/bin/bash

i=1
for PARAM in srr dl_path output_path
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done

mkdir -p $output_path
cd $output_path
wget $dl_path/$srr;
/bioinfo/local/build/sratoolkit.2.2.2b/bin/fastq-dump -I --split-files $srr;
sed -E "s/^((@|\+)SRR[^.]+\.[^.]+)\.(1|2)/\1/" $srr\_1.fastq > $srr\_1.fixed.fastq;
sed -E "s/^((@|\+)SRR[^.]+\.[^.]+)\.(1|2)/\1/" $srr\_2.fastq > $srr\_2.fixed.fastq;