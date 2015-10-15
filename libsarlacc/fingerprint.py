# Core imports
import os
import sys

# Library imports
import numpy as np

# Local imports
from .calc import (
        spearman_roc,
        kendall_tau,
        absolute_distance,
        get_dist_mat,
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

metrics = {'sp': spearman_roc,
           'kt': kendall_tau,
           'hd': absolute_distance}

def process_files(files, png=False, metric='sp', output=None):
    resolution = 100
    reader = DataFileReader({'d_e':'d_e', 'd_i':'d_i'},
                            FingerprintData)

    descriptors = batch_process(files, reader)

    if len(descriptors) < 2:
        log("Need at least 2 things to compare!", cat='error')
        return

    histograms, names = zip(*[(x.histogram, x.name) for x in descriptors])
    mat = get_dist_mat(histograms, metric=metrics[metric])

    clusters = cluster(mat, names,
                            method='centroid',
                            distance=0.4)
    logClosestPair(mat, names)
    logFarthestPair(mat, names)

    if output:
        write_mat_file(output,
                       mat,
                       np.array(names, dtype='S10'),
                       clusters)
    if png:
        prefix = os.path.curdir
        log('Writing files to {}'.format(prefix))
        #plot_file(x, y, fname=outfile, nbins=resolution)

