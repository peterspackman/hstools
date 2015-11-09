# Core imports
import os

# Library imports
import numpy as np
import pandas as pd
from hdbscan import HDBSCAN

# Local imports
from .calc import cluster
from .config import log, logClosestPair, logFarthestPair
from .datafile import (
        batch_process,
        write_mat_file,
        DataFileReader,
        HarmonicsData
        )

shape_keys = {'radius':'radius', 'coefficients':'coefficients', 'invariants':'invariants'}
dnorm_keys = {'radius':'radius', 'coefficients':'dnorm_coefficients', 'invariants':'dnorm_invariants'}
curvature_keys = {'radius':'radius', 'coefficients':'curvature_coefficients', 'invariants':'curvature_invariants'}

modes = {'shape':shape_keys, 'dnorm':dnorm_keys, 'curvature':curvature_keys}

def process_files(files, no_radius=False, mode='shape', output=None, **kwargs):
    dendrogram = None
    method = HDBSCAN
    distance = 0.4
    use_radius = True
    log('Reading {} files...'.format(len(files)))

    reader = DataFileReader(modes[mode], HarmonicsData)

    descriptors = batch_process(files, reader, procs=1)

    if len(descriptors) < 2:
        log("Need at least 2 things to compare!", cat='error')
        return

    radius, invariants, names  = zip(*[(x.radius, x.invariants, x.name.stem) for x in descriptors])
    columns = [x for x in range(0,len(invariants[0]))]
    if not no_radius:
        columns = ['r'] + columns
        invariants = [np.append(r, x) for r, x in zip(radius, invariants)]
    mat = np.array(invariants)
    df = pd.DataFrame(mat, columns=columns)
    clusters = cluster(mat, method=method, min_cluster_size=5)
    df['name'] = names
    df['cluster'] = clusters
    if output:
        write_mat_file(output, mat, names, clusters)
    #log('{}'.format(df[['name', 'cluster']]))

