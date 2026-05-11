#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=1:mem=5gb
# Change the path of `short_username` as needed
#PBS -o /rds/general/user/{short_username}/home/EPIONCHO.IBM/logs/
#PBS -e /rds/general/user/{short_username}/home/EPIONCHO.IBM/logs/

filetorun="${FILETORUN}" # change this as needed
# outputname="${}" # change this as needed

export OUTPUT_PREFIX="${OUTPUTPREFIX}"
export OUTPUT_FOLDER="${OUTPUTFOLDERNAME}"

echo "${OUTPUT_PREFIX}"
echo "${OUTPUT_FOLDER}"

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source ~/anaconda3/etc/profile.d/conda.sh

cp "$HOME/EPIONCHO.IBM/${filetorun}" $TMPDIR
mkdir -p $TMPDIR/inst/extdata
mkdir $TMPDIR/rfils
mkdir $TMPDIR/rout
mkdir $TMPDIR/$OUTPUT_FOLDER

cp $HOME/EPIONCHO.IBM/inst/extdata/eye_disease_probabilties_updated.rds $TMPDIR/inst/extdata

mkdir -p $HOME/EPIONCHO.IBM/rfils
mkdir -p $HOME/EPIONCHO.IBM/rout

# Change to the $TMPDIR
cd $TMPDIR

echo "R is about to run"

R CMD BATCH "$TMPDIR/${filetorun}" "$TMPDIR/rout/${filetorun}Output_${PBS_ARRAY_INDEX}.txt"

echo "Copying R runtime files to ${HOME}/EPIONCHO.IBM/rout"
cp -r $TMPDIR/rfils $HOME/EPIONCHO.IBM/
cp -r $TMPDIR/rout $HOME/EPIONCHO.IBM/

echo "Copying R output files to ${HOME}/EPIONCHO.IBM/output/"
cp $TMPDIR/${OUTPUT_FOLDER}/${OUTPUTPREFIX}_${PBS_ARRAY_INDEX}.rds $HOME/EPIONCHO.IBM/output/


echo "R has finished running"