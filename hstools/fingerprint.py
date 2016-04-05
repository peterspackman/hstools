# Core imports
import os
import sys

# Library imports
import numpy as np
from tqdm import tqdm

# Local imports
from .calc import cluster
from .config import log
from .datafile import (
        batch_process,
        write_mat_file,
        DataFileReader,
        FingerprintData,
        hexbin_plotfile
)


def flatten_hist(h):
    """
    Helper method for extracting the flattened
    array of a 2D histogram
    """
    return h[0].flatten()


def process_files(files, png=False, metric='sp', output=None):
    """
    Given a list of hdf5 files (Path objects), create the Hirshfeld
    fingerprint from the data, then perform a clustering on the results.

    Allows Hirshfeld fingerprints to be saved as part of the process.

    Returns df, a pandas.DataFrame object containing the data
    """

    resolution = 5
    reader = DataFileReader({'d_e': 'd_e', 'd_i': 'd_i'},
                            FingerprintData)

    descriptors = batch_process(files, reader)

    if len(descriptors) < 2:
        log("Need at least 2 things to compare!", cat='error')
        return

    histograms, names = zip(*[(x.histogram, x.name) for x in descriptors])
    mat = np.array(list((map(flatten_hist, histograms))))
    clusters = cluster(mat)
    if output:
        write_mat_file(output,
                       mat,
                       np.array(names, dtype='S10'),
                       clusters)

    if png:
        prefix = os.path.curdir
        log('Writing fingerprint plot files to {}'.format(prefix))
        for c in tqdm(descriptors, unit='figure'):
            outfile = str(c.name) + '.png'
            hexbin_plotfile(c.d_i, c.d_e, fname=outfile)
