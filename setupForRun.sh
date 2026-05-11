#!/bin/bash
folder="EPIONCHO.IBM"
outputfolder="output"
file=""
helpval="Usage: $0 [-f file to run] [-n number of runs] [-c clean output and logs folder] [-m model folder (default: EPIONCHO.IBM)] [-o output folder (default: output)] [-p output prefix (defailt: test_)]"
clean=false
numberofruns=200
outputprefix="test_"
while getopts ":f:m:o:n:c:p:" opt; do
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

if [[ -z "$file" ]]; then
  echo "Error: Option -f is required." >&2
  echo $helpval >&2
  exit 1
fi

if [[ "$clean" == true ]]; then
    echo "TODO: Should be cleaning"
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

numrunoption="1-${numberofruns}"
if [[ $numberofruns == *-* ]]; then
  numrunoption="${numberofruns}"
fi

echo "qsub -J value: ${numrunoption}"
echo "output prefix: ${outputprefix}"
echo "output folder name: ${outputfolder}"

if [ $(pwd) == "$HOME/${folder}" ]; then
  qsub -J "$numrunoption" -v "OUTPUTFOLDERNAME=${outputfolder},FILETORUN=${file},OUTPUTPREFIX=${outputprefix}" run.sh
else
   echo "Not in the right folder"
fi
