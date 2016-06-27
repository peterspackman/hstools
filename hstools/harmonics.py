"""
Command line functionality for clustering HS based on
spherical harmonic shape descriptors
"""
# Library imports
import numpy as np
import pandas as pd
from hdbscan import HDBSCAN

# Local imports
from .calc import cluster
from .config import log
from .datafile import (
    batch_process,
    DataFileReader,
    HarmonicsData
)

SHAPE_KEYS = {'radius': 'radius',
              'coefficients': 'coefficients',
              'invariants': 'invariants'}

DNORM_KEYS = {'radius': 'radius',
              'coefficients': 'dnorm_coefficients',
              'invariants': 'dnorm_invariants'}

MODES = {'shape': SHAPE_KEYS, 'dnorm': DNORM_KEYS}

def make_dataframe(invariants, names, clusters, columns):
    """ Construct a dataframe for the invariants, names, clusters """
    output_data = pd.DataFrame(invariants, columns=columns)
    output_data['name'] = names
    output_data['cluster'] = clusters

def process_files(files, no_radius=False, mode='shape', **kwargs):
    """
    Given a list of hdf5 files (Path objects), read the spherical harmonics
    shape descriptors data, then perform a clustering on the results.

    Returns output_data, a pandas.DataFrame object containing the data
    """
    #dendrogram = kwargs.get('dendrogram', None)
    method = kwargs.get('method', HDBSCAN)
    #distance = kwargs.get('distance', 0.4)

    log('Reading {} files...'.format(len(files)))

    reader = DataFileReader(MODES[mode], HarmonicsData)

    descriptors = batch_process(files, reader, procs=1)

    if len(descriptors) < 2:
        log("Need at least 2 things to compare!", cat='error')
        return

    radius, invariants, names = zip(*[(x.radius,
                                       x.invariants,
                                       x.name) for x in descriptors])

    columns = [x for x in range(0, len(invariants[0]))]
    if not no_radius:
        columns = ['r'] + columns
        invariants = [np.append(r, x) for r, x in zip(radius, invariants)]
    invariants = np.array(invariants)
    clusters = cluster(invariants, method=method, min_cluster_size=5)
    return make_dataframe(invariants, names, clusters, columns)
