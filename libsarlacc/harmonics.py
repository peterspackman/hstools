# Core imports
import os

# Library imports
import numpy as np

# Local imports
from .calc import get_dist_mat, euclidean, cluster
from .config import log, logClosestPair, logFarthestPair
from .datafile import (
        batch_process,
        write_mat_file,
        DataFileReader,
        HarmonicsData
        )

shape_keys = {'coefficients':'coefficients', 'invariants':'invariants'}
dnorm_keys = {'coefficients':'dnorm_coefficients', 'invariants':'dnorm_invariants'}
curvature_keys = {'coefficients':'curvature_coefficients', 'invariants':'curvature_invariants'}

modes = {'shape':shape_keys, 'dnorm':dnorm_keys, 'curvature':curvature_keys}

def process_files(files, mode='shape', output=None):
    metric = euclidean
    dendrogram = None
    method = 'ward'
    distance = 0.4
    output = None

    reader = DataFileReader(modes[mode], HarmonicsData)

    descriptors = batch_process(files, reader, procs=1)

    if len(descriptors) < 2:
        log("Need at least 2 things to compare!", cat='error')
        return

    invariants, names  = zip(*[(x.invariants, x.name) for x in descriptors])

    mat = get_dist_mat(invariants, metric=metric)
    clusters = cluster(mat, names, dendrogram=dendrogram,
                            method=method, distance=distance)

    if output:
        write_mat_file(output,
                       mat,
                       np.array(names, dtype='S10'),
                       clusters)

    logClosestPair(mat, names)
    logFarthestPair(mat, names)
