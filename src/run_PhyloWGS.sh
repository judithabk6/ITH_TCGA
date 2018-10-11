i=1
for PARAM in patient folder pythonenv phyloWGS_path smchet_path
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done


source $pythonenv/bin/activate

rm -rf results/$patient/PhyloWGS/$folder
./src/format_cnv_phylowgs.py $patient


mkdir -p results/$patient/PhyloWGS/$folder
python src/create_phylowgs_inputs.py --vcf-type S1=csv_pyclone --cnvs S1=results/$patient/PhyloWGS/cnv_table.csv --output-cnvs results/$patient/PhyloWGS/$folder/parsed_cnv_data.txt --output-variants results/$patient/PhyloWGS/$folder/parsed_ssm_data.txt --output-params results/$patient/PhyloWGS/$folder/params.json S1=results/$patient/pyclone/$folder/input.tsv
cd results/$patient/PhyloWGS/$folder
python $phyloWGS_path/evolve.py parsed_ssm_data.txt parsed_cnv_data.txt


python $phyloWGS_path/write_results.py --include-ssm-names $patient\_$folder trees.zip summ.json.gz muts.json.gz mutass.zip

mkdir outputs
export PYTHONPATH=$phyloWGS_path:$PYTHONPATH
python $smchet_path/smchet-challenge/create-smchet-report/write_report.py summ.json.gz muts.json.gz mutass.zip outputs
cd ../../../..

