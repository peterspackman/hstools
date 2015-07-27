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
from docopt import docopt
import numpy as np
import scipy.sparse

# Local imports
from . import calc
from .data import log, log_traceback, log_error,  \
    logClosestPair, logFarthestPair
from .fileio import proc_file_hist, batch_hist, write_mat_file

test_f = {'sp': calc.spearman_roc,
          'kt': calc.kendall_tau,
          'hd': calc.absolute_distance}

test_names = {'sp': 'Spearman ROC',
              'kt': "Kendall's Tau",
              'hd': 'absolute histogram distance'}


def test_from_arg(key):
    if key not in test_f:
        log_error('{} is not a valid option.'.format(key))
        sys.exit(1)
    return test_f[key]

def process_file_list(files, args, procs):
    dendrogram = args['--dendrogram']
    method = args['--method']
    distance = float(args['--distance'])
    test = test_from_arg(args['--test'])

    if len(files) < 1:
        log_error('No files to process.')
        sys.exit(1)

    descriptors = batch_hist(files,
                             resolution=int(args['--bins']),
                             save_figs=args['--save-figures'],
                             procs=procs)

    if len(descriptors) < 2:
        log_error("Need at least 2 things to compare!")
        return

    histograms, names = zip(*[(x.histogram, x.name) for x in descriptors])
    mat = calc.get_dist_mat(histograms, test=test,
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
    mtest = calc.euclidean
    files = []

    for file_pattern in args['<filepattern>']:
        if os.path.isfile(file_pattern):
            fname = file_pattern
            files.append(file_pattern)

        elif os.path.isdir(file_pattern):
            from .fileio import glob_directory
            new_files = glob_directory(file_pattern, '*.h5')
            files += new_files
            if len(new_files) < 1:
                log_error("{} no files matching {}.".format(file_pattern, '*.h5'))
        else:
            from glob import glob
            new_files = glob(file_pattern)
            files += new_files
            if len(new_files) < 1:
                log_error("pattern {} matched no files.".format(file_pattern))

    if len(files) > 0:
        process_file_list(files, args, procs)
    else:
        log_error('No files found...')
