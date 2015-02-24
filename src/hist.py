"""
Usage:
    sarlacc hist [options] [<filepattern>...]

    -h, --help
    -n, --dry-run
    -v, --verbose
    -b=NUM, --bins=NUM             Set the number of bins to use
                                   in the histogram. [default: 100]
    -t=TEST, --test=TEST           Select which test will be used.
                                   [default: sp]
    -p, --save-figures             Plot histograms calculated and
                                   save them to file.
    -o=FILE, --output=FILE         Write the result to a given file.
                                   Works for both surface and hist modes.
                                   Will save the distance matrix in hist mode,
                                   and write S.A. infor in surface mode
    -d=FILE, --dendrogram=FILE     Save the the generated clustering as a
                                   dendrogram.
    -m=METHOD, --method=METHOD     Use METHOD when calculating linkage. One of
                                   'average', 'single', 'complete', 'weighted',
                                   'centroid', 'median', 'ward'.
                                   [default: complete]
    --distance=THRESHOLD           The threshold distance for leaves to be
                                   classified as clustered. Unlikely to change
                                   much. [default: 0.4]
"""

# Core imports
import os
import re
import sys
# Library imports
import numpy as np
import scipy.sparse
from docopt import docopt
from .data import log, log_traceback, log_error, logClosestPair, logFarthestPair
from . import calc
from .fileio import proc_file_hist, batch_hist, write_mat_file


test_f = {'sp': calc.spearman_roc,
          'kt': calc.kendall_tau,
          'hd': calc.absolute_distance}

test_names = {'sp': 'Spearman rank order coefficient',
              'kt': "Kendall's Tau",
              'hd': 'Sigma histogram distance'}


def process_file_list(files, args, procs):
    dendrogram = args['--dendrogram']
    method = args['--method']
    distance = float(args['--distance'])
    if len(files) < 1:
        log_error('No files to process.')
        sys.exit(1)

    histograms, names = batch_hist(files,
                                   resolution=int(args['--bins']),
                                   save_figs=args['--save-figures'],
                                   procs=procs)

    mat = calc.get_dist_mat(histograms, test=test_f[args['--test']],
                            threads=procs*2)

    clusters = calc.cluster(mat, names,
                            dendrogram=dendrogram,
                            method=method,
                            distance=distance)
    logClosestPair(mat, names)
    logFarthestPair(mat, names)

    if args['--output']:
        fname = args['--output']
        write_mat_file(fname,
                       mat,
                       np.array(names, dtype='S10'),
                       clusters)


def hist_main(argv, procs=4):
    args = docopt(__doc__, argv=argv)

    if len(args['<filepattern>']) < 2:
        file_pattern = args['<filepattern>'][0]

        if os.path.isfile(file_pattern):
            fname = file_pattern
            if not save_figs:
                log('Not saving figure, so this '
                    'command will have no output')
            h, name = proc_file_hist(fname, resolution=args['--bins'],
                                     save_figs=args['--save_figs'])

        elif os.path.isdir(file_pattern):
            from .fileio import glob_directory
            with glob_directory(file_pattern, '*.hdf5') as files:
                process_file_list(files, args, procs)

    else:
        process_file_list(args['<filepattern>'], args, procs)
