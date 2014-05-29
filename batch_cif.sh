#!/bin/sh
# This script is designed to run multiple instances of tonto
# in order to calculate .cxs files from .cif inputs

LOG=batch_cif.log
DIR=$1
parallel  --joblog $LOG --progress '~/phd/run_cif_to_surface.exe -res 0.5 -cif {} -cxs {.}.cxs > {.}.out' ::: $DIR/*.cif
