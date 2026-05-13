#!/bin/bash
folder="EPIONCHO.IBM"
outputfolder="output"
file="all_funcs_combined.R"
postprocessfile="process_multiple_runs.R"
helpval="Usage: $0 [-f file to run (default: all_funcs_combined.R)] [-n number of runs] [-c clean output and logs folder] [-m model folder (default: EPIONCHO.IBM)] [-o output folder (default: output)] [-p output prefix (defailt: test)] [-g path to save processed data. Must include filename, and end in .rds] [-b run model and post process]"
clean=false
numberofruns=200
outputprefix="test"
processdatapath=""
runandpostprocess=false
while getopts ":f:m:o:n:cp:g:b" opt; do
  case $opt in
    b) 
      echo "Running model and post-processing"
      runandpostprocess=true
      ;;
    c) 
      echo "Cleaning folders before running (this may take a long time)"
      clean=true
      ;;
    f) file="$OPTARG" ;;
    n) numberofruns="$OPTARG" ;;
    m) folder="$OPTARG" ;;
    o) outputfolder="$OPTARG" ;;
    p) outputprefix="$OPTARG" ;;
    g) processdatapath="$OPTARG" ;;
    :)
      if [[ "$OPTARG" == "f" ]]; then
        echo "Error: Option -f requires an argument." >&2
        echo $helpval >&2
        exit 1
      fi
      ;;
    \?)
      echo $helpval >&2
      exit 1
      ;;
  esac
done

if [[ -n "${processdatapath}" ]]; then
  if [[ "$runandpostprocess" == false ]]; then
    numberofruns="1"
  fi
  mkdir -p $(dirname "${processdatapath}")
  if ! [ -d "${outputfolder}" ] && [ "$runandpostprocess" == false ]; then
    echo "Folder containing outputs does not exist."
    exit 1
  fi
fi

if [[ -z "$file" ]]; then
  echo "Error: Option -f is required." >&2
  echo $helpval >&2
  exit 1
fi

if [[ "$clean" == true ]]; then
    echo "Cleaning logs/"
    rm -rf "$HOME/${folder}/logs/"
    echo "Cleaning rout/$file"
    rm -rf "$HOME/${folder}/rout/$file"*
    exit 0
fi

if ! [ -d "${outputfolder}" ]; then
    echo "No existing output folder, making a new one"
    mkdir -p $HOME/${folder}/${outputfolder}/
fi
if ! [ -d "logs" ]; then
    echo "No existing logs folder, making a new one"
    mkdir -p $HOME/${folder}/logs/
fi

if ! [ -d "rout" ]; then
    echo "No existing rout folder, making a new one"
    mkdir -p $HOME/${folder}/rout/
fi

if [ -f "$HOME/miniforge3/bin/conda" ]; then
  echo "conda is installed in miniforge3: $(~/miniforge3/bin/conda --version)"
else
  echo "conda not found in miniforge3, loading miniforge3 and setting up conda"
  module load miniforge/3
  miniforge-setup
fi

eval "$(~/miniforge3/bin/conda shell.bash hook)"
source ~/miniforge3/etc/profile.d/conda.sh

if conda env list | grep -q "^my_r_env "; then
    echo "environment 'my_r_env' exists"
else
    echo "environment 'my_r_env' not found, creating it now."
    # R version likely could be updated, this is the last version the 
    # model was verified on
    conda create -n my_r_env r-base=4.4.3 -c conda-forge
fi
conda activate my_r_env 

dos2unix $HOME/EPIONCHO.IBM/${file}

dos2unix $HOME/EPIONCHO.IBM/run.sh

echo "output prefix: ${outputprefix}"
echo "output folder name: ${outputfolder}"

numrunoption="1-${numberofruns}"
if [[ $numberofruns == *-* ]]; then
  numrunoption="${numberofruns}"
fi

if [ $(pwd) == "$HOME/${folder}" ]; then
  if [[ -n "${processdatapath}" && "$runandpostprocess" == false ]]; then
    echo "qsub for processing outputs"
    job_id=$(qsub -v "OUTPUTFOLDERNAME=${outputfolder},FILETORUN=${postprocessfile},OUTPUTPREFIX=${outputprefix},MODELFOLDER=${folder},PROCESSDATAPATH=${processdatapath}" \
      -o $HOME/$folder/logs/ -e $HOME/$folder/logs/ \
      run_process_files.sh)
  elif [[ $numberofruns == "1" ]]; then
    echo "qsub for single run"
    job_id=$(qsub -v "OUTPUTFOLDERNAME=${outputfolder},FILETORUN=${file},OUTPUTPREFIX=${outputprefix},MODELFOLDER=${folder}" \
      -o $HOME/$folder/logs/ -e $HOME/$folder/logs/ \
      run.sh)
    exit 0
  else
    job_id=$(qsub -J "$numrunoption" -v "OUTPUTFOLDERNAME=${outputfolder},FILETORUN=${file},OUTPUTPREFIX=${outputprefix},MODELFOLDER=${folder}" \
      -o $HOME/$folder/logs/ -e $HOME/$folder/logs/ \
      run.sh)
  fi
  echo "Job ID: $job_id"
  if [[ "$runandpostprocess" == true ]]; then 
    dashWstring="depend=afterok:$job_id"
    echo "qsub for processing outputs"
    new_job_id=$(qsub -v "OUTPUTFOLDERNAME=${outputfolder},FILETORUN=${postprocessfile},OUTPUTPREFIX=${outputprefix},MODELFOLDER=${folder},PROCESSDATAPATH=${processdatapath}" \
      -o $HOME/$folder/logs/ -e $HOME/$folder/logs/ \
      -W $dashWstring \
      run_process_files.sh)
    echo "Post Process Job ID: $new_job_id"
  fi
else
   echo "Not in the right folder"
fi
