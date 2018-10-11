#cd $PBS_O_WORKDIR

i=1
for PARAM in patient folder Rpath
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i | sed 's/ /,/g') ;
done


export PATH=$Rpath:$PATH
./src/run_sciclone.R -p $patient -f $folder