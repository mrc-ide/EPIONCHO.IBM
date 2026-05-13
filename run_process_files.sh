#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=1:mem=60gb


filetorun="${FILETORUN}"
export OUTPUT_PREFIX="${OUTPUTPREFIX}"
export OUTPUT_FOLDER="${OUTPUTFOLDERNAME}"
export PROCESS_DATA_PATH="${PROCESSDATAPATH}"
echo "${OUTPUT_PREFIX}"
echo "${OUTPUT_FOLDER}"
echo "${PROCESS_DATA_PATH}"

eval "$(~/miniforge3/bin/conda shell.bash hook)"
source ~/miniforge3/etc/profile.d/conda.sh
conda activate my_r_env

mkdir -p $HOME/${MODELFOLDER}/rfils
mkdir -p $HOME/${MODELFOLDER}/rout

cd $HOME/${MODELFOLDER}
echo "R is about to run"

R CMD BATCH "${filetorun}" "$HOME/${MODELFOLDER}/rout/${filetorun}.ROutput.txt"

echo "R has finished running"