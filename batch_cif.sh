#!/bin/sh
LOG=batch_cif.log
DIR=$1


parallel  --joblog $LOG --progress '~/phd/run_cif_to_surface.exe -cif {} -cxs {.}.cxs > {.}.out' ::: $DIR/*.cif
