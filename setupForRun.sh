#!/bin/bash
folder="EPIONCHO.IBM"
outputfolder="output"
file="all_funcs_combined.R"
helpval="Usage: $0 [-f file to run (default: all_funcs_combined.R)] [-n number of runs] [-c clean output and logs folder] [-m model folder (default: EPIONCHO.IBM)] [-o output folder (default: output)] [-p output prefix (defailt: test_)] [-g folder to save processed data]"
clean=false
numberofruns=200
outputprefix="test"
processdatafolder=""
while getopts ":f:m:o:n:c:p:g:" opt; do
  case $opt in
    c) 
      echo "Cleaning folders before running"
      clean=true
      ;;
    f) file="$OPTARG" ;;
    n) numberofruns="$OPTARG" ;;
    m) folder="$OPTARG" ;;
    o) outputfolder="$OPTARG" ;;
    p) outputprefix="$OPTARG" ;;
    g) processdatafolder="$OPTARG" ;;
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

if [[ -n "${processdatafolder}" ]]; then
  file="process_multiple_runs.R"
  numberofruns="1"
  mkdir -p ${processdatafolder}
  if ! [ -d "${outputfolder}" ]; then
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
    rm -rf "$HOME/${folder}/logs/$file"*
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

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source ~/anaconda3/etc/profile.d/conda.sh

conda activate myenv 

dos2unix $HOME/EPIONCHO.IBM/${file}

dos2unix $HOME/EPIONCHO.IBM/run.sh

echo "output prefix: ${outputprefix}"
echo "output folder name: ${outputfolder}"

numrunoption="1-${numberofruns}"
if [[ $numberofruns == *-* ]]; then
  numrunoption="${numberofruns}"
fi

if [ $(pwd) == "$HOME/${folder}" ]; then
  if [[ -n "${processdatafolder}" ]]; then
    echo "qsub for processing outputs"
    qsub -v "OUTPUTFOLDERNAME=${outputfolder},FILETORUN=${file},OUTPUTPREFIX=${outputprefix},MODELFOLDER=${folder},PROCESSDATAFOLDER=${processdatafolder}" \
      -o $HOME/$folder/logs/ -e $HOME/$folder/logs/ \
      run_process_files.sh
    exit 0
  fi
  if [[ $numberofruns == "1" ]]; then
    echo "qsub for single run"
    qsub -v "OUTPUTFOLDERNAME=${outputfolder},FILETORUN=${file},OUTPUTPREFIX=${outputprefix},MODELFOLDER=${folder},PROCESSDATAFOLDER=${processdatafolder}" \
      -o $HOME/$MODELFOLDER/logs/ -e $HOME/$folder/logs/ \
      run.sh
    exit 0
  fi
  qsub -J "$numrunoption" -v "OUTPUTFOLDERNAME=${outputfolder},FILETORUN=${file},OUTPUTPREFIX=${outputprefix},MODELFOLDER=${folder}" \
    -o $HOME/$folder/logs/ -e $HOME/$folder/logs/ \
    run.sh
  exit 0
else
   echo "Not in the right folder"
fi
