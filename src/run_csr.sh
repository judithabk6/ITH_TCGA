#!/bin/bash

i=1
for PARAM in patient
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done

source config.txt
source /data/users/jabecass/dl_tools_centos/miniconda3_activate
source activate my_python36_centos
export PATH=$RPATH:$PATH

./src/write_csr_input.py $patient
cwd=$(pwd)

for folder in protected_hg38_vcf public_hg38_vcf protected_hg38_vcf_absolute;  do
rm -rf results/$patient/CSR/$folder/input
rm -rf results/$patient/CSR/$folder/output
rm -rf results/$patient/CSR/$folder/preprocessed

mkdir -p results/$patient/CSR/$folder/input
mkdir -p results/$patient/CSR/$folder/output
mkdir -p results/$patient/CSR/$folder/preprocessed
run_ok=true
for ith_method in PhyloWGS pyclone sciclone expands
do 
if [ -f results/$patient/$ith_method/$folder/input_csr.tsv ]; then ln -s $cwd/results/$patient/$ith_method/$folder/input_csr.tsv $cwd/results/$patient/CSR/$folder/input/$ith_method\_input_csr.tsv;
else run_ok=false ; echo $run_ok; echo $ith_method; fi done

if [ "$run_ok" = true ] ; then
    Rscript src/CSR_preprocess.R $cwd/results/$patient/CSR/$folder/input $cwd/results/$patient/CSR/$folder/preprocessed 5000 0.1
    python src/CSR.py $cwd/results/$patient/CSR/$folder/preprocessed $cwd/results/$patient/CSR/$folder/output -500
fi
done


