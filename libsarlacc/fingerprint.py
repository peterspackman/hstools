# Core imports
import os
import sys

# Library imports
import numpy as np

# Local imports
from .calc import (
        cluster,
)

from .config import (
        log,
        logClosestPair,
        logFarthestPair
)

from .datafile import (
        batch_process,
        write_mat_file,
        DataFileReader,
        FingerprintData
)


def process_files(files, png=False, metric='sp', output=None):
    resolution = 5
    reader = DataFileReader({'d_e':'d_e', 'd_i':'d_i'},
                            FingerprintData)

    descriptors = batch_process(files, reader)

    if len(descriptors) < 2:
        log("Need at least 2 things to compare!", cat='error')
        return

    histograms, names = zip(*[(x.histogram, x.name.stem) for x in descriptors])
    flatten_hist = lambda h: h[0].flatten()
    mat = np.array(list((map(flatten_hist, histograms))))
    print(mat)

    clusters = cluster(mat) 
    if output:
        write_mat_file(output,
                       mat,
                       np.array(names, dtype='S10'),
                       clusters)
    if png:
        prefix = os.path.curdir
        log('Writing files to {}'.format(prefix))
        #plot_file(x, y, fname=outfile, nbins=resolution)

