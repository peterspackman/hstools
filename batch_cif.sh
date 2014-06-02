#!/bin/bash
# This script is designed to run multiple instances of tonto
# in order to calculate .cxs files from .cif inputs
EXE=/home/prs/phd/run_cif_to_surface.exe
LOG=batch_cif.log
RES=0.5
OPTIND=1

while getopts "h?l:r:" opt; do
  case "$opt" in
    h|\?)
      echo "help"
      exit 0
      ;;
    r)
      RES=$OPTARG
      ;;
    l)
      LOG=$OPTARG
      ;;
  esac
done
shift $(($OPTIND - 1))
DIR=$1
echo "Processing .cxs files in $DIR with resolution of $RES"

### DO THE JOBS ###

JOB="$EXE -res $RES -cif {} -cxs {.}.cxs > {.}.out"
parallel  --joblog $LOG --progress $JOB ::: $DIR/*.cif
