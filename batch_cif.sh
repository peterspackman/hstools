#!/bin/sh
# This script is designed to run multiple instances of tonto
# in order to calculate .cxs files from .cif inputs
EXE=/home/prs/phd/run_cif_to_surface.exe
LOG=batch_cif.log
RES=0.25
DIR=$1
JOB="$EXE -res $RES -cif {} -cxs {.}.cxs > {.}.out"

parallel  --joblog $LOG --progress $JOB ::: $DIR/*.cif
