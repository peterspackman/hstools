#!/bin/sh
DIR=cif
parallel $* --progress '~/phd/run_cif_to_surface.exe -cif {} -cxs {.}.cxs > {.}.out' ::: $DIR/*.cif
