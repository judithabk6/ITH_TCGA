#!/usr/bin/env bash

source config.txt 

mkdir $DATA_FOLDER
cd $DATA_FOLDER
mkdir BRCA BLCA HNSC
for loc in BRCA BLCA HNSC; do
cd $loc;
mkdir $CLINICAL_FOLDERNAME  $PROTECTED_FOLDERNAME  $PUBLIC_FOLDERNAME;
cd ..;
done
cd ..


mkdir tmp_dl
cd tmp_dl
mkdir clinical_data
cd clinical_data
wget http://download.cbioportal.org/brca_tcga.tar.gz
wget http://download.cbioportal.org/blca_tcga.tar.gz
wget http://download.cbioportal.org/hnsc_tcga.tar.gz
mkdir BRCA BLCA HNSC
cd BLCA
tar -xvf ../blca_tcga.tar.gz
cd ..
cd BRCA
tar -xvf ../brca_tcga.tar.gz
cd ..
cd HNSC
tar -xvf ../hnsc_tcga.tar.gz
cd ..
for loc in BRCA BLCA HNSC; do
ln -s $PWD/$loc/data_bcr_clinical_data_patient.txt ../../$DATA_FOLDER/$loc/$CLINICAL_FOLDERNAME/data_bcr_clinical_data_patient.txt ;
ln -s $PWD/$loc/data_bcr_clinical_data_sample.txt ../../$DATA_FOLDER/$loc/$CLINICAL_FOLDERNAME/data_bcr_clinical_data_sample.txt ;
done
cd ..
mkdir mutation_data
cd mutation_data
../../gdc-client download -m ../../gdc_manifest.2018-07-23.txt -t ../../$GDC_TOKEN
for folder in `ls .`; do 
if [[ $folder = *".txt" ]]; then
echo "just the manifest"
else
mutation_file=`ls $folder/*.gz`
#mutation_file=`ls e$folder/*.maf`
cancer_loc=`echo $mutation_file | cut -d "/" -f 2 | cut -d "." -f 2`;
gunzip $mutation_file
full_filename=`ls $folder/*.maf`
filename=`echo $full_filename | cut -d "/" -f 2`;
protection_status=`echo $mutation_file | cut -d "/" -f 2 | cut -d "." -f 7`;
if [ "$protection_status" = "protected" ]; then 
    ln -s $PWD/$full_filename ../../$DATA_FOLDER/$cancer_loc/$PROTECTED_FOLDERNAME/$filename
elif [ "$protection_status" = "somatic" ]; then 
    ln -s $PWD/$full_filename ../../$DATA_FOLDER/$cancer_loc/$PUBLIC_FOLDERNAME/$filename
fi
if [ ! -f "../../$DATA_FOLDER/$cancer_loc/annotations.txt" ] ; then
    ln -s $PWD/$folder/annotations.txt ../../$DATA_FOLDER/$cancer_loc/annotations.txt
fi
fi
done



