#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=1:mem=60gb


filetorun="${FILETORUN}"
export OUTPUT_PREFIX="${OUTPUTPREFIX}"
export OUTPUT_FOLDER="${OUTPUTFOLDERNAME}"
export PROCESS_DATA_FOLDER="${PROCESSDATAFOLDER}"
echo "${OUTPUT_PREFIX}"
echo "${OUTPUT_FOLDER}"
echo "${PROCESS_DATA_FOLDER}"


eval "$(~/anaconda3/bin/conda shell.bash hook)"
source ~/anaconda3/etc/profile.d/conda.sh
conda activate myenv 

mkdir -p $HOME/${MODELFOLDER}/rfils
mkdir -p $HOME/${MODELFOLDER}/rout

cd $HOME/${MODELFOLDER}
echo "R is about to run"

R CMD BATCH "${filetorun}" "$HOME/${MODELFOLDER}/rout/${filetorun}.ROutput.txt"

echo "R has finished running"