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
    values, names = batch_harmonics(files, procs=procs, property=args['--property'])
    if len(values) < 2:
        log_error("Need at least 2 things to compare!")
        return
    coefficients, invariants = zip(*values)
    mat = calc.get_dist_mat(invariants, test=mtest, threads=procs*2)
    clusters = calc.cluster(mat, names, dendrogram=dendrogram,
                            method=method, distance=distance)

    if args['--output']:
        fname = args['--output']
        write_mat_file(fname,
                       mat,
                       np.array(names, dtype='S10'),
                       clusters)

    logClosestPair(mat, names)
    logFarthestPair(mat, names)


def harmonics_main(argv, procs=4):
    args = docopt(__doc__, argv=argv)
    mtest = calc.euclidean

    if len(args['<filepattern>']) < 2:

        file_pattern = args['<filepattern>'][0]

        if os.path.isfile(file_pattern):
            fname = file_pattern
            values, cname = proc_file_harmonics(fname, property=args['--property'])
            coefficients, invariants = values
            log(cname)
            log(coefficients)

        elif os.path.isdir(file_pattern):
            from .fileio import glob_directory
            with glob_directory(file_pattern, '*.hdf5') as files:
                process_file_list(files, args, procs)
        else:
            from glob import glob
            files = glob(file_pattern)
            if len(files) < 1:
                log_error("pattern {} matched no files.".format(file_pattern))
            else:
                process_file_list(files, args, procs)

    else:
        process_file_list(args['<filepattern>'], args, procs)
