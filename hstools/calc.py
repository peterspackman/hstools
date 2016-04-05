"""
Collection of calculations based on HS data
"""

# Core imports
from collections import OrderedDict

# Library imports
import numpy as np
import scipy.stats
from hdbscan import HDBSCAN
from periodictable import elements

# Local imports
from .config import log


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
    """Return the interaction at [i,j] as a string e.g. C->H
    """
    symbols = (str(elements[i+1]), str(elements[j+1]))
    if order:
        return "{0} -> {1}".format(*symbols)
    else:
        return "{0} -> {1}".format(*sorted(symbols))


def get_contrib(surface_area, order=False):
    """ Given a the triangles that make up a hirshfeld surface,
    and lists of the closest internal and external atoms along
    with their respective distances from the surface,
    calculate the makeup of the hirshfeld surface in terms of
    which element->element interactions are responsible for that
    area """
    contributions = OrderedDict()
    contributions_percent = OrderedDict()
    # setting defaults for these

    for i, j in np.transpose(np.nonzero(surface_area)):
        # Key in the form "internal -> external" e.g. "F -> H"
        key = interaction(i, j, order=order)
        if key not in contributions:
            contributions[key] = 0.0
        contributions[key] += surface_area[i, j]

    for contrib in contributions:
        contributions_percent[contrib] = \
                np.round(contributions[contrib] / sum(contributions.values()),
                         decimals=8)
    return contributions, contributions_percent


def bin_data(xvals, yvals, bins=10, bounds=False):
    """ Puts the data x & y into a given number of bins.
        Currently, bounds simply tells the program to use
        set bins."""

    if not bounds:
        xmin, xmax = np.min(xvals), np.max(xvals)
        ymin, ymax = np.min(yvals), np.max(yvals)
    else:
        xmin, ymin = 0.5, 0.5
        xmax, ymax = 2.5, 2.5

    diffx = (xmax - xmin) / (bins - 1.0)
    diffy = (ymax - ymin) / (bins - 1.0)

    weights = np.ones(len(xvals))

    # this is a slightly modified version of np.digitize()
    xyi = np.vstack((xvals, yvals)).T
    xyi -= [xmin, ymin]
    xyi /= [diffx, diffy]
    xyi = np.floor(xyi, xyi).T

    grid = scipy.sparse.coo_matrix((weights, xyi), shape=(bins, bins)).toarray()

    return grid, np.linspace(xmin, xmax, bins), np.linspace(ymin, ymax, bins)
