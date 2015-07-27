"""
Usage:
    sarlacc harmonics [options] [<filepattern>...]

    -h, --help
    -n, --dry-run
    -v, --verbose
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
    -p=NAME, --property=NAME       Use invariants of property on surface. On of
                                   'dnorm', 'curvature', 'shape'.
                                   [default: shape]
"""
# Core imports
import os

# Library imports
from docopt import docopt
import numpy as np

# Local imports
from . import calc
from .data import log, log_error, logClosestPair, logFarthestPair
from .fileio import proc_file_harmonics, batch_harmonics, write_mat_file


def process_file_list(files, args, procs):
    mtest = calc.euclidean
    dendrogram = args['--dendrogram']
    method = args['--method']
    distance = float(args['--distance'])
    descriptors = batch_harmonics(files,
                                  procs=procs,
                                  property=args['--property'])
    if len(descriptors) < 2:
        log_error("Need at least 2 things to compare!")
        return

    invariants, names  = zip(*[(x.invariants, x.name) for x in descriptors])

    mat = calc.get_dist_mat(invariants, test=mtest, threads=procs*2)
    clusters = calc.cluster(mat, names, dendrogram=dendrogram,
                            method=method, distance=distance)

    if args['--output']:
        write_mat_file(args['--output'],
                       mat,
                       np.array(names, dtype='S10'),
                       clusters)

    logClosestPair(mat, names)
    logFarthestPair(mat, names)


def harmonics_main(argv, procs=4):
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
