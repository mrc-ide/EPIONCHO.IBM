#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=5gb

filetorun="${FILETORUN}" # change this as needed
# outputname="${}" # change this as needed

export OUTPUT_PREFIX="${OUTPUTPREFIX}"
export OUTPUT_FOLDER="${OUTPUTFOLDERNAME}"

echo "${OUTPUT_PREFIX}"
echo "${OUTPUT_FOLDER}"

eval "$(~/miniforge3/bin/conda shell.bash hook)"
source ~/miniforge3/etc/profile.d/conda.sh
conda activate my_r_env 

cp "$HOME/${MODELFOLDER}/${filetorun}" $TMPDIR
mkdir -p $TMPDIR/inst/extdata
mkdir $TMPDIR/rfils
mkdir $TMPDIR/rout
mkdir $TMPDIR/$OUTPUT_FOLDER

cp $HOME/${MODELFOLDER}/inst/extdata/eye_disease_probabilties_updated.rds $TMPDIR/inst/extdata

mkdir -p $HOME/${MODELFOLDER}/rfils
mkdir -p $HOME/${MODELFOLDER}/rout

# Change to the $TMPDIR
cd $TMPDIR

echo "R is about to run"

R CMD BATCH "$TMPDIR/${filetorun}" "$TMPDIR/rout/${filetorun}Output_${PBS_ARRAY_INDEX}.txt"

echo "Copying R runtime files to ${HOME}/${MODELFOLDER}/rout"
cp -r $TMPDIR/rfils $HOME/${MODELFOLDER}/
cp -r $TMPDIR/rout $HOME/${MODELFOLDER}/

echo "Copying R output files to ${HOME}/${MODELFOLDER}/${OUTPUT_FOLDER}/"
cp $TMPDIR/${OUTPUT_FOLDER}/${OUTPUT_PREFIX}_${PBS_ARRAY_INDEX}.rds $HOME/${MODELFOLDER}/${OUTPUT_FOLDER}


echo "R has finished running"