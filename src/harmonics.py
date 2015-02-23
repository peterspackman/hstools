"""
Usage:
    sarlacc harmonics [options] [<filepattern>]

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
"""
import os

from docopt import docopt
import numpy as np
from . import calc
from .data import log
from .fileio import proc_file_harmonics, batch_harmonics, write_mat_file
from .modes import logClosestPair, logFarthestPair


def harmonics_main(argv, procs=4):
    args = docopt(__doc__, argv=argv)
    mtest = calc.euclidean

    if os.path.isfile(args['<filepattern>']):
        fname = args['<filepattern>']
        values, cname = proc_file_harmonics(fname)
        coefficients, invariants = values
        log(cname)
        log(coefficients)

    elif os.path.isdir(args['<filepattern>']):
        dendrogram = args['--dendrogram']
        method = args['--method']
        distance = float(args['--distance'])
        dirname = args['<filepattern>']
        values, names = batch_harmonics(dirname, procs=procs)
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
