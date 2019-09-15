#!/bin/bash

source config.txt 
###########################################
### step 1. Download and organize files ###
###########################################
mkdir logs
source $python35_env_path/bin/activate
./src/organize_data_download.sh
echo 1,BRCA >> 20180801_get_RNAseq.csv
echo 2,BLCA >> 20180801_get_RNAseq.csv
echo 3,HNSC >> 20180801_get_RNAseq.csv

qsub -t 1-3 -N 20180801_get_RNAseq -q batch -d $PWD -l walltime=15:00:00,mem=30gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20180801_get_RNAseq.csv src/get_RNAseq_data.sh
deactivate
# in python 3.6 environnement
source $python36_env_path/bin/activate
# ASCAT copy number
segment_liftover -l $liftover_path -i ASCAT_TCGA -o liftover_ASCAT_TCGA -c hg19ToHg38 -si BLCA.ascat.segments.txt -so BLCA.ascat.segments.hg38.txt
segment_liftover -l $liftover_path -i ASCAT_TCGA -o liftover_ASCAT_TCGA -c hg19ToHg38 -si BRCA.ascat.segments.txt -so BRCA.ascat.segments.hg38.txt
segment_liftover -l $liftover_path -i ASCAT_TCGA -o liftover_ASCAT_TCGA -c hg19ToHg38 -si HNSC.ascat.segments.txt -so HNSC.ascat.segments.hg38.txt

# ABSOLUTE copy number
./src/get_input_absolute_liftover.py
segment_liftover -l $liftover_path -i external_data -o external_data -c hg19ToHg38 -si TCGA_mastercalls.abs_segtabs.fixed.liftover_input.txt -so TCGA_mastercalls.abs_segtabs.fixed.liftover_output.txt

./src/get_cn_formatted_liftover.py
deactivate

rm -f 20180807_get_snv_filters.csv
j=1
for cancer_loc in BRCA BLCA HNSC ; do
echo $j,$cancer_loc,$python35_env_path >> 20180807_get_snv_filters.csv
j=$((j+1)); 
done

qsub -t 1-3 -N 20180807_get_snv_filters -q batch -d $PWD -l walltime=10:00:00,mem=60gb,nodes=1:ppn=1 -v INPUT_FILE=20180807_get_snv_filters.csv src/get_SNV_filters_and_data_exploration.sh
./src/get_final_snv_files_merge.py
mv 20180807_get_snv_filters.e* 20180807_get_snv_filters.o* logs

################################
### step 2. Run ITH methods  ###
################################
# get pyclone input
cat official_patient_list.csv |  while read a; do  echo "$a,public_hg38_vcf,tmp/useful_final_public_merge_cnv_purity.csv" >> 20180801_public_input.csv;   done
qsub -t 1-1790 -N 20180801_public_input_pyclone -q batch -d $PWD -l walltime=0:20:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20180801_public_input.csv src/build_pyclone_input.sh
mv 20180801_public_input_pyclone.e* 20180801_public_input_pyclone.o* logs

cat official_patient_list.csv |  while read a; do  echo "$a,protected_hg38_vcf,tmp/useful_final_qc_merge_cnv_purity.csv" >> 20180801_protected_input.csv;   done
qsub -t 1-1790 -N 20180801_protected_input_pyclone -q batch -d $PWD -l walltime=0:20:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20180801_protected_input.csv src/build_pyclone_input.sh
mv 20180801_protected_input_pyclone.e* 20180801_protected_input_pyclone.o* logs

cat official_patient_list.csv |  while read a; do  echo "$a,protected_hg38_vcf_absolute,tmp/useful_final_qc_merge_cnv_purity_absolute.csv" >> 20190801_protected_absolute_input.csv;   done
qsub -t 1-1790%98 -N 20190801_protected_absolute_input_pyclone -q batch -d $PWD -l walltime=0:20:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20190801_protected_absolute_input.csv src/build_pyclone_input.sh
mv 20190801_protected_absolute_input_pyclone.e* 20190801_protected_absolute_input_pyclone.o* logs

# pyclone
cat official_patient_list.csv |  while read a; do  echo "$a,public_hg38_vcf,$python27_env_path" >> 20180801_public_run_pyclone.csv;   done
qsub -t 1-1790 -N 20180801_public_run_pyclone -q batch -d $PWD -l walltime=0:20:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20180801_public_run_pyclone.csv src/run_pyclone.sh

cat official_patient_list.csv |  while read a; do  echo "$a,protected_hg38_vcf,$python27_env_path" >> 20180801_protected_run_pyclone.csv;   done
qsub -t 1-1790 -N 20180801_protected_run_pyclone -q batch -d $PWD -l walltime=0:20:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20180801_protected_run_pyclone.csv src/run_pyclone.sh
mv 20180801_public_run_pyclone.e* 20180801_public_run_pyclone.o* 20180801_protected_run_pyclone.o* 20180801_protected_run_pyclone.e* logs


cat official_patient_list.csv |  while read a; do  echo "$a,protected_hg38_vcf_absolute,$python27_env_path" >> 20190801_protected_run_pyclone.csv;   done
qsub -t 1-1790%98 -N 20190801_protected_run_pyclone -q batch -d $PWD -l walltime=0:20:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20190801_protected_run_pyclone.csv src/run_pyclone.sh
mv 20190801_protected_run_pyclone.o* 20190801_protected_run_pyclone.e* logs



# sciclone
cat official_patient_list.csv |  while read a; do  echo "$a,protected_hg38_vcf,$RPATH" >> 20180801_protected_hg38_sciclone.csv;   done
cat official_patient_list.csv |  while read a; do  echo "$a,protected_hg38_vcf_absolute,$RPATH" >> 20190801_protected_absolute_hg38_sciclone.csv;   done
cat official_patient_list.csv |  while read a; do  echo "$a,public_hg38_vcf,$RPATH" >> 20180801_public_hg38_sciclone.csv;   done
qsub -t 1-1790 -N 20180801_protected_hg38_sciclone -q batch -d $PWD -l walltime=1:00:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20180801_protected_hg38_sciclone.csv src/run_sciclone.sh
qsub -t 1-1790 -N 20180801_public_hg38_sciclone -q batch -d $PWD -l walltime=1:00:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20180801_public_hg38_sciclone.csv src/run_sciclone.sh
mv 20180801_protected_hg38_sciclone.e* 20180801_protected_hg38_sciclone.o* 20180801_public_hg38_sciclone.e* 20180801_public_hg38_sciclone.o* logs
qsub -t 1-1790%90 -N 20190801_protected_absolute_hg38_sciclone -q batch -d $PWD -l walltime=1:00:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20190801_protected_absolute_hg38_sciclone.csv src/run_sciclone.sh
mv 20190801_protected_absolute_hg38_sciclone.e* 20190801_protected_absolute_hg38_sciclone.o* logs

# baseline - same input file as pyclone input
qsub -t 1-1790 -N 20180801_baseline_protected -q batch -d $PWD -l walltime=0:20:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20180801_protected_input.csv src/run_baseline.sh
qsub -t 1-1790 -N 20180801_baseline_public -q batch -d $PWD -l walltime=0:20:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20180801_public_input.csv src/run_baseline.sh
mv 20180801_baseline_protected.e* 20180801_baseline_protected.o* 20180801_baseline_public.e* 20180801_baseline_public.o* logs

# expands
qsub -t 1-1790 -N 20180801_run_expands -q batch -d $PWD -l walltime=10:00:00,mem=1gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20180801_public_hg38_sciclone.csv src/run_expands.sh
qsub -t 1-1790 -N 20180801_run_expands_protected -q batch -d $PWD -l walltime=10:00:00,mem=1gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20180801_protected_hg38_sciclone.csv src/run_expands.sh
qsub -t 1-1790%92 -N 20190801_run_expands_protected_absolute -q batch -d $PWD -l walltime=10:00:00,mem=1gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20190801_protected_absolute_hg38_sciclone.csv src/run_expands.sh
mv 20180801_run_expands.e* 20180801_run_expands.o* 20180801_run_expands_protected.e* 20180801_run_expands_protected.o* logs
mv 20190801_run_expands_protected_absolute.e* 20190801_run_expands_protected_absolute.o* logs

# phyloWGS
cat official_patient_list.csv |  while read a; do     echo "$a,protected_hg38_vcf,$python27_env_path,$phyloWGS_path,$smchet_path" >> 20180801_run_phylowgs_protected.csv;    done
cat official_patient_list.csv |  while read a; do     echo "$a,protected_hg38_vcf_absolute,$python27_env_path,$phyloWGS_path,$smchet_path" >> 20190801_run_phylowgs_protected_absolute.csv;    done
cat official_patient_list.csv |  while read a; do     echo "$a,public_hg38_vcf,$python27_env_path,$phyloWGS_path,$smchet_path" >> 20180801_run_phylowgs_public.csv;    done
qsub -t 1-1790 -N 20180801_run_phylowgs_public -q batch -d $PWD -l walltime=10:00:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20180801_run_phylowgs_public.csv src/run_PhyloWGS.sh
qsub -t 1-1790 -N 20180801_run_phylowgs_protected -q batch -d $PWD -l walltime=10:00:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20180801_run_phylowgs_protected.csv src/run_PhyloWGS.sh
qsub -t 1-5 -N 20190801_run_phylowgs_protected_absolute -q batch -d $PWD -l walltime=10:00:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20190801_run_phylowgs_protected_absolute.csv src/run_PhyloWGS.sh

mv 20180801_run_phylowgs_public.e* 20180801_run_phylowgs_public.o* 20180801_run_phylowgs_protected.e* 20180801_run_phylowgs_protected.o* logs

qsub -t 1-1790 -N 20180801_csr -q batch -d $PWD -l walltime=0:20:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=official_patient_list.csv src/run_csr.sh
mv 20180801_csr.e* 20180801_csr.o* logs

source $python35_env_path/bin/activate
./src/get_feature_table.py $TORQUE_LOG_PATH
./src/perform_survival_regression.py
./src/perform_survcomp_stats.R
./src/get_article_graphics.py
deactivate
