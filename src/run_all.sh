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
qsub -t 6-1790%90 -N 20190801_run_phylowgs_protected_absolute -q batch -d $PWD -l walltime=10:00:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20190801_run_phylowgs_protected_absolute.csv src/run_PhyloWGS.sh
mv 20180801_run_phylowgs_public.e* 20180801_run_phylowgs_public.o* 20180801_run_phylowgs_protected.e* 20180801_run_phylowgs_protected.o* logs

# run CSR
qsub -t 1-1790%80 -N 20190801_csr -q batch -d $PWD -l walltime=0:30:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=official_patient_list.csv src/run_csr.sh
mv 20190801_csr.e* 20190801_csr.o* logs

# get additional metrics on samples for which all methods ran
qsub -t 1-690%180 -N 20190912_metrics_ITH -q batch -d $PWD -l walltime=00:20:00,mem=3gb,nodes=1:ppn=1 -v INPUT_FILE=intersection_patient_list.csv src/get_standard_metrics_ith.sh

##################################
### step 3. get final results  ###
##################################
source $python35_env_path/bin/activate
./src/get_feature_table.py $TORQUE_LOG_PATH
./src/perform_survival_regression.py
./src/perform_survcomp_stats.R
./src/get_article_graphics.py
deactivate

# actually the get_feature_table was this time run on the cluster
qsub -N get_feature_table -q batch -d $PWD -l walltime=10:00:00,mem=10gb,nodes=1:ppn=1 src/get_feature_table.sh


#####################################
### step 4. single cell analysis  ###
#####################################
# download single cell datasets 
# (first from SRA, and then for some datasets from the ENA platform of the EBI
# directly in fastq format)
qsub -t 1-14 -N 20190829_dl_sra -q batch -d $PWD -l walltime=0:05:00,mem=3gb,nodes=1:ppn=1  -v INPUT_FILE=20190829_dl_single_cell_dataset.csv src/download_srr.sh
qsub -t 1-16 -N 20190906_dl_sra -q batch -d $PWD -l walltime=10:00:00,mem=3gb,nodes=1:ppn=1  -v INPUT_FILE=20190906_dl_crc_embl.csv src/download_crc.sh
qsub -N 20190829_tnbc_fastq -l walltime=20:00:00,mem=3gb,nodes=1:ppn=1 download_tnbc_fastq.sh 


# run the WES pipeline
# the code is not mine to share
# TNBC run
ASSEMBLY=hg19
OUTPUT_CONFIG=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/single_cell_dataset/BRCA_wang/ewok_config.txt
TARGET_BED=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/external_data/SeqCap_EZ_Exome_v3_hg19_capture_targets.bed

/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/ewok-devel/ewok.sh GetConf --template  ewok-devel/CONFIG_TEMPLATE --configFile ewok-devel/species_design_configs.csv --designType wes --genomeAssembly ${ASSEMBLY} --outputConfig ${OUTPUT_CONFIG}  --targetBed ${TARGET_BED}

NAME=BRCA_wang
OUTPUT_DIR=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/single_cell_dataset/BRCA_wang/ewok_results
OUTPUT_CONFIG=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/single_cell_dataset/BRCA_wang/ewok_config.txt
SAMPLE_PLAN=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/20190902_BRCA_wang_TNBC_ewok_plan.csv
RUN=BRCA_wang
/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/ewok-devel/ewok_run.sh  --conf ${OUTPUT_CONFIG}  --samplePlan ${SAMPLE_PLAN} --outDir ${OUTPUT_DIR}  --run ${RUN}

# CRC CO5
ASSEMBLY=hg19
OUTPUT_CONFIG=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/single_cell_dataset/CRC_leung/ewok_config_CO5.txt
TARGET_BED=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/external_data/TruSeq_exome_targeted_regions.hg19.bed
/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/ewok-devel/ewok.sh GetConf --template  ewok-devel/CONFIG_TEMPLATE --configFile ewok-devel/species_design_configs.csv --designType wes --genomeAssembly ${ASSEMBLY} --outputConfig ${OUTPUT_CONFIG}  --targetBed ${TARGET_BED}

SAMPLE_PLAN=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/20190903_CRC_CO5_sample_plan.csv
OUTPUT_DIR=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/single_cell_dataset/CRC_leung/ewok_results
RUN=CRC_leung_CO5

/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/ewok-devel/ewok_run.sh  --conf ${OUTPUT_CONFIG}  --samplePlan ${SAMPLE_PLAN} --outDir ${OUTPUT_DIR}  --run ${RUN}

# CRC CO8
ASSEMBLY=hg19
OUTPUT_CONFIG=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/single_cell_dataset/CRC_leung/ewok_config_CO8.txt
TARGET_BED=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/external_data/SeqCap_EZ_Exome_v3_hg19_capture_targets.bed
/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/ewok-devel/ewok.sh GetConf --template  ewok-devel/CONFIG_TEMPLATE --configFile ewok-devel/species_design_configs.csv --designType wes --genomeAssembly ${ASSEMBLY} --outputConfig ${OUTPUT_CONFIG}  --targetBed ${TARGET_BED}

SAMPLE_PLAN=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/20190903_CRC_sample_plan.csv
OUTPUT_DIR=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/single_cell_dataset/CRC_leung/ewok_results
RUN=CRC_leung_CO8

/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/ewok-devel/ewok_run.sh  --conf ${OUTPUT_CONFIG}  --samplePlan ${SAMPLE_PLAN} --outDir ${OUTPUT_DIR}  --run ${RUN}

# ALL data
ASSEMBLY=hg19
OUTPUT_CONFIG=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/single_cell_dataset/ALL_gawad/ewok_config.txt
TARGET_BED=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/external_data/nexterarapidcapture_expandedexome_uniqueintervals.bed
TARGET_BED=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/external_data/all_genome.bed

/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/ewok-devel/ewok.sh GetConf --template  ewok-devel/CONFIG_TEMPLATE --configFile ewok-devel/species_design_configs.csv --designType wes --genomeAssembly ${ASSEMBLY} --outputConfig ${OUTPUT_CONFIG}  --targetBed ${TARGET_BED}

SAMPLE_PLAN=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/20190903_ALL_sample_plan.csv
OUTPUT_DIR=/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/single_cell_dataset/ALL_gawad/ewok_results
RUN=ALL_gawad
/data/users/jabecass/ITHeterogene/rerun_tcga_het_project/ewok-devel/ewok_run.sh  --conf ${OUTPUT_CONFIG}  --samplePlan ${SAMPLE_PLAN} --outDir ${OUTPUT_DIR}  --run ${RUN}


# run DepthOfCoverage to get complementary mutations
JAVA=/bioinfo/local/build/Centos/java/jdk1.8.0_74/bin;
GATKV=/bioinfo/local/build/Centos/GATK/GenomeAnalysisTK-3.5;
$JAVA/java -Xmx10g -jar $GATKV/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R /data/annotations/Galaxy/Human/hg19/bwa/hg19.fasta \
-I single_cell_dataset/ALL_gawad/ewok_results/ALL_gawad/SRR1517761/Preprocess/SRR1517761_uniq.nodup.onTarget.q20.recal.bam \
-o single_cell_dataset/ALL_gawad/doc_SRR1517761.doc \
-baseCounts \
-L external_data/gawad_ground_truth_p1.list \
--omitIntervalStatistics --includeDeletions --minMappingQuality 20 --minBaseQuality 13

$JAVA/java -Xmx10g -jar $GATKV/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R /data/annotations/Galaxy/Human/hg19/bwa/hg19.fasta \
-I single_cell_dataset/ALL_gawad/ewok_results/ALL_gawad/SRR1517762/Preprocess/SRR1517762_uniq.nodup.onTarget.q20.recal.bam \
-o single_cell_dataset/ALL_gawad/doc_SRR1517762.doc \
-baseCounts \
-L external_data/gawad_ground_truth_p1.list \
--omitIntervalStatistics --includeDeletions --minMappingQuality 20 --minBaseQuality 13

$JAVA/java -Xmx10g -jar $GATKV/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R /data/annotations/Galaxy/Human/hg19/bwa/hg19.fasta \
-I single_cell_dataset/ALL_gawad/ewok_results/ALL_gawad/SRR1517763/Preprocess/SRR1517763_uniq.nodup.onTarget.q20.recal.bam \
-o single_cell_dataset/ALL_gawad/doc_SRR1517763.doc \
-baseCounts \
-L external_data/gawad_ground_truth_p2.list \
--omitIntervalStatistics --includeDeletions --minMappingQuality 20 --minBaseQuality 13

$JAVA/java -Xmx10g -jar $GATKV/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R /data/annotations/Galaxy/Human/hg19/bwa/hg19.fasta \
-I single_cell_dataset/ALL_gawad/ewok_results/ALL_gawad/SRR1517764/Preprocess/SRR1517764_uniq.nodup.onTarget.q20.recal.bam \
-o single_cell_dataset/ALL_gawad/doc_SRR1517764.doc \
-baseCounts \
-L external_data/gawad_ground_truth_p2.list \
--omitIntervalStatistics --includeDeletions --minMappingQuality 20 --minBaseQuality 13


for srr in SRR3472566 SRR3472567 SRR3472569 SRR3472571; do
JAVA=/bioinfo/local/build/Centos/java/jdk1.8.0_74/bin;
GATKV=/bioinfo/local/build/Centos/GATK/GenomeAnalysisTK-3.5;
$JAVA/java -Xmx10g -jar $GATKV/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R /data/annotations/Galaxy/Human/hg19/bwa/hg19.fasta \
-I single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO5/$srr/Preprocess/$srr\_uniq.nodup.onTarget.q20.recal.bam \
-o single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO5/doc_$srr.doc \
-baseCounts \
-L external_data/CO5_malilik_positions.list \
--omitIntervalStatistics --includeDeletions --minMappingQuality 20 --minBaseQuality 13
done

for srr in SRR3472796 SRR3472798 SRR3472799 SRR3472800; do
JAVA=/bioinfo/local/build/Centos/java/jdk1.8.0_74/bin;
GATKV=/bioinfo/local/build/Centos/GATK/GenomeAnalysisTK-3.5;
$JAVA/java -Xmx10g -jar $GATKV/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R /data/annotations/Galaxy/Human/hg19/bwa/hg19.fasta \
-I single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO8/$srr/Preprocess/$srr\_uniq.nodup.onTarget.q20.recal.bam \
-o single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO8/doc_$srr.doc \
-baseCounts \
-L external_data/CO8_malilik_positions.list \
--omitIntervalStatistics --includeDeletions --minMappingQuality 20 --minBaseQuality 13
done

for srr in SRR1163508 SRR1298936; do
JAVA=/bioinfo/local/build/Centos/java/jdk1.8.0_74/bin;
GATKV=/bioinfo/local/build/Centos/GATK/GenomeAnalysisTK-3.5;
$JAVA/java -Xmx10g -jar $GATKV/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R /data/annotations/Galaxy/Human/hg19/bwa/hg19.fasta \
-I single_cell_dataset/BRCA_wang/ewok_results/BRCA_wang/$srr/Preprocess/$srr\_uniq.nodup.onTarget.q20.recal.bam \
-o single_cell_dataset/BRCA_wang/ewok_results/BRCA_wang/doc_$srr.doc \
-baseCounts \
-L external_data/tnbc_positions_liftover_malilik.list \
--omitIntervalStatistics --includeDeletions --minMappingQuality 20 --minBaseQuality 13
done

./src/get_single_cell_ith_method_inpu.py

# run ITH methods on the single cell dataset
rm 20180801_singlecell_protected_run_pyclone.csv
cat single_cell_patient_list.csv |  while read a; do  echo "$a,protected_hg38_vcf,$python27_env_path" >> 20180801_singlecell_protected_run_pyclone.csv;   done
qsub -t 1-7 -N 20180801_singlecell_protected_run_pyclone -q batch -d $PWD -l walltime=02:00:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20180801_singlecell_protected_run_pyclone.csv src/run_pyclone.sh

rm 20180801_singlecell_protected_hg38_sciclone.csv
cat single_cell_patient_list.csv |  while read a; do  echo "$a,protected_hg38_vcf,$RPATH" >> 20180801_singlecell_protected_hg38_sciclone.csv;   done
qsub -t 1-7 -N 20180801_singlecell_protected_hg38_sciclone -q batch -d $PWD -l walltime=1:00:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20180801_singlecell_protected_hg38_sciclone.csv src/run_sciclone.sh
qsub -t 1-7 -N 20180801_singlecell_run_expands_protected -q batch -d $PWD -l walltime=10:00:00,mem=1gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20180801_singlecell_protected_hg38_sciclone.csv src/run_expands.sh

rm 20180801_singlecell_run_phylowgs_protected.csv
cat single_cell_patient_list.csv |  while read a; do     echo "$a,protected_hg38_vcf,$python27_env_path,$phyloWGS_path,$smchet_path" >> 20180801_singlecell_run_phylowgs_protected.csv;    done
qsub -t 1-7 -N 20180801_singlecell_run_phylowgs_protected -q batch -d $PWD -l walltime=10:00:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=20180801_singlecell_run_phylowgs_protected.csv src/run_PhyloWGS.sh

qsub -t 1-7 -N 20190801_csr_single_cell -q batch -d $PWD -l walltime=0:30:00,mem=3gb,nodes=1:ppn=1 -m ae -M judithabk6@gmail.com -v INPUT_FILE=single_cell_patient_list.csv src/run_csr.sh

./src/get_single_cell_metrics.py
