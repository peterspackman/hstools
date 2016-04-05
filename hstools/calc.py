#!/usr/bin/python

# Core imports
from collections import defaultdict

# Library imports
from matplotlib import pyplot as plt
import numpy as np
import scipy.stats
from hdbscan import HDBSCAN
from periodictable import elements

# Local imports
from .config import log, Timer


def cluster(raw_data, method=HDBSCAN, **kwargs):
    """
    Takes a 2D array of raw data (observations)
    and performs the desired clustering on the
    data.

    By default will use HDBSCAN.

    Returns clustering information.

    """
    log('Clustering using {}'.format(method.__name__))
    clusterer = method(**kwargs)
    return clusterer.fit_predict(raw_data)


def interaction(i, j, order=False):
    symbols = (str(elements[i+1]), str(elements[j+1]))
    if order:
        return "{0} -> {1}".format(*symbols)
    else:
        return "{0} -> {1}".format(*sorted(symbols))


def get_contrib(sa, order=False):
    """ Given a the triangles that make up a hirshfeld surface,
    and lists of the closest internal and external atoms along
    with their respective distances from the surface,
    calculate the makeup of the hirshfeld surface in terms of
    which element->element interactions are responsible for that
    area """
    contrib = defaultdict(float)
    contrib_p = defaultdict(float)
    # setting defaults for these

    for i, j in np.transpose(np.nonzero(sa)):
        # Key in the form "internal -> external" e.g. "F -> H"
        key = interaction(i, j, order=order)
        contrib[key] += sa[i, j]

    for x in contrib:
        p = np.round(contrib[x] / sum(contrib.values()), decimals=8)
        contrib_p[x] = p

    return contrib, contrib_p


def bin_data(x, y, bins=10, bounds=False):
    """ Puts the data x & y into a given number of bins.
        Currently, bounds simply tells the program to use
        set bins."""

    nx = ny = bins

    if not bounds:
        xmin, xmax = min(x), max(x)
        ymin, ymax = min(y), max(y)
    else:
        xmin, ymin = 0.5, 0.5
        xmax, ymax = 2.5, 2.5

    dx = (xmax - xmin) / (nx - 1.0)
    dy = (ymax - ymin) / (ny - 1.0)

    weights = np.ones(len(x))

    # this is a slightly modified version of np.digitize()
    xyi = np.vstack((x, y)).T
    xyi -= [xmin, ymin]
    xyi /= [dx, dy]
    xyi = np.floor(xyi, xyi).T

    grid = scipy.sparse.coo_matrix((weights, xyi), shape=(nx, ny)).toarray()

    return grid, np.linspace(xmin, xmax, nx), np.linspace(ymin, ymax, ny)
