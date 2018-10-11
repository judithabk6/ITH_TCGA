

i=1
for PARAM in loc
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i | sed 's/ /,/g') ;
done

source config.txt
export PATH=$RPATH33:$PATH
mkdir -p data/$loc/RNAseq
./src/get_RNAseq_data.R $loc
