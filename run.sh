#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=1:mem=5gb
# Change the path of `short_username` as needed
#PBS -o /rds/general/user/{short_username}/home/EPIONCHO.IBM/logs/
#PBS -e /rds/general/user/{short_username}/home/EPIONCHO.IBM/logs/

filetorun="run_model_oae.R" # change this as needed
outputname="oem_test" # change this as needed

module load intel-suite
module load anaconda3/personal

cp "$HOME/EPIONCHO.IBM/${filetorun}" $TMPDIR
mkdir $TMPDIR/data
mkdir $TMPDIR/rfils
mkdir $TMPDIR/rout
mkdir $TMPDIR/output
cp $HOME/EPIONCHO.IBM/data/eye_disease_probabilties_updated.rds $TMPDIR/data

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
cp $TMPDIR/output/${outputname}_${PBS_ARRAY_INDEX}.rds $HOME/EPIONCHO.IBM/output/


echo "R has finished running"