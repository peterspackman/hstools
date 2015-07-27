"""
Usage:
    sarlacc corr [options] [<filepattern>...]

    -h, --help
    -n, --dry-run
    -v, --verbose
    -o=FILE, --output=FILE         Write the result to a given file.
                                   Works for both surface and hist modes.
                                   Will save the distance matrix in hist mode,
                                   and write S.A. infor in surface mode
"""
import os
# Library imports
from docopt import docopt
import numpy as np

from .fileio import readh5file
from .data import log


def print_correlation(filename):
    data = readh5file(filename, ['vertices', 'd_norm', 'indices'])
    vertices, dnorm, faces = data['vertices'], data['d_norm'], data['indices']
    center = np.apply_along_axis(np.mean, 0, vertices)
    print("Centered at: {}".format(center))
    get_radius = lambda x: np.linalg.norm(x - center)
    radii = np.apply_along_axis(get_radius, 1, vertices)
    radii = radii/np.mean(radii)
    log(np.corrcoef(radii, dnorm))


def main(argv, procs=4):
    args = docopt(__doc__, argv=argv)

    if len(args['<filepattern>']) < 2:

        file_pattern = args['<filepattern>'][0]

        if os.path.isfile(file_pattern):
            fname = file_pattern
            print_correlation(fname)

        elif os.path.isdir(file_pattern):
            from .fileio import glob_directory
            with glob_directory(file_pattern, '*.h5') as files:
                for f in files:
                    print_correlation(f)
        else:
            from glob import glob
            files = glob(file_pattern)
            if len(files) < 1:
                log_error("pattern {} matched no files.".format(file_pattern))
            else:
                for f in files:
                    print_correlation(f)
