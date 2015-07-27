"""
Usage:
    sarlacc clusters [options] [<filepattern>...]

    -h, --help
    -n, --dry-run
    -v, --verbose
    -o=FILE, --output=FILE         Write the result to a given file.
                                   Works for both surface and hist modes.
                                   Will save the distance matrix in hist mode,
                                   and write S.A. infor in surface mode
    -m=METHOD, --method=METHOD     Only use METHOD when checking. One of
                                   'average', 'single', 'complete', 'weighted',
                                   'centroid', 'median', 'ward'.
                                   [default: complete]
    --distance=THRESHOLD           The threshold distance for leaves to be
                                   classified as clustered. Unlikely to change
                                   much. [default: 0.4]
    -p=NAME, --property=NAME       Only calculate for property NAME. On of
                                   'dnorm', 'curvature', 'shape'.
                                   [default: shape]
"""
# Core imports
import os

# Library imports
from docopt import docopt
import numpy as np

# Local imports
from .harmonics import process_file_list as harmonics_list
from .hist import process_file_list as hist_list
from . import calc
from .data import log, log_error, logClosestPair, logFarthestPair
from .fileio import get_files_from_pattern


methods = ['average', 'single', 'complete', 'weighted',
           'centroid', 'median', 'ward']

properties = ['curvature', 'shape', 'dnorm']


def check_clusters(files):
    log(get_files_from_pattern(files))


def clusters_main(argv, procs=4):
    args = docopt(__doc__, argv=argv)

    if len(args['<filepattern>']) < 2:

        file_pattern = args['<filepattern>'][0]

        if os.path.isfile(file_pattern):
            log_error('Please supply a directory or list of files.')
            sys.exit(1)

        elif os.path.isdir(file_pattern):
            check_clusters(file_pattern)
        else:
            if len(files) < 1:
                log_error("pattern {} matched no files.".format(file_pattern))
            else:
                check_clusters(file_pattern)
    else:
        check_clusters(args['<filepattern>'])
