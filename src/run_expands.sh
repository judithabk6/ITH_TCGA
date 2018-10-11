#cd $PBS_O_WORKDIR

i=1
for PARAM in patient folder Rpath
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i | sed 's/ /,/g') ;
done

source config.txt
export JAVA_HOME="$java_path/"
export JAVA_JRE="$java_path/jre/bin"
export PATH="$java_path/bin:$PATH"
export PATH="$java_path/jre/bin:$PATH"
export JAVA_CPPFLAGS="-I$java_path/include -I$java_path/include/linux"
export JAVA="$JAVA_HOME/bin/java"
export JAVA_LIBS="-L$java_path/jre/lib/amd64/server -ljvm"
export LD_LIBRARY_PATH=$JAVA_HOME/jre/lib/amd64:$JAVA_HOME/jre/lib/amd64/server


export PATH=$Rpath:$PATH
./src/run_expands.R -p $patient -f $folder