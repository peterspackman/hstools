#!/usr/bin/python
# Core imports
from collections import defaultdict
from functools import reduce
from itertools import combinations
import concurrent.futures
import json
# Library imports
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram as dend
import fastcluster as fc
import numpy as np
import progressbar as pb
import scipy.cluster.hierarchy
import scipy.spatial.distance
import scipy.stats as stats
# Local imports
from . import data
from .data import log, logger, Timer, elements


def spearman_roc(histograms):
    """ Calculate the Spearman rank-order correlation coefficient from
    2 histograms This may need to be modified, I'm uncertain whether or
    not 2 zeroes are ignored READ likely that they aren't, as such
    artificially high correlations are probable"""
    H1, H2 = histograms
    hist1, _, _ = H1
    hist2, _, _ = H2
    x = hist1.flatten()
    y = hist2.flatten()

    r, p = stats.spearmanr(x, y)
    return r


def kendall_tau(histograms):
    """ Calculate Kendall's Tau from the given histograms. Significantly slower
    than Spearman ROC, and seems to produce slightly worse results."""
    H1, H2 = histograms
    hist1, _, _ = H1
    hist2, _, _ = H2
    x = hist1.flatten()
    y = hist2.flatten()
    r, p = stats.kendalltau(x, y)
    return r


def absolute_distance(histograms):
    """ Calculate the absolute distance between two histograms"""
    H1, H2 = histograms
    x, _, _ = H1
    y, _, _ = H2
    d = np.sum(np.subtract(x, y))
    return abs(d)


def euclidean(x):
    i1, i2 = x
    d = np.power(np.sum(np.power(np.subtract(i2, i1), 2)), 0.5)
    return d


# insert a vector of length l into a new matrix,
# keeping values from mat into the upper triangle
def insert_into_distance_mat(mat, vec):
    tmp = np.zeros((vec.size, vec.size))
    tmp[0, :] = vec
    tmp[1:, 1:] = mat
    return tmp


def write_dendrogram_file(fname, Z, names, no_labels=False,
                          test_name='test', distance=0.4):
    if not distance:
        distance = 0.4
    threshold = distance * np.max(Z[:, 2])
    # Plot stuff
    plt.xlabel('Compound')
    plt.ylabel('Dissimilarity')
    plt.style.use('ggplot')
    dpi = 200
    plt.suptitle("""Clustering dendrogram of {0}
                compounds using {1}""".format(len(names), test_name))
        # Create a dendrogram
    dend(Z, no_labels=no_labels, color_threshold=threshold, labels=names)
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=90)
    if len(names) > 100:
        fig = plt.gcf()
        fig.set_size_inches(10.5, min(len(names)*0.1, 32768/dpi))
    plt.savefig(fname, dpi=dpi, bbox_inches='tight')
    plt.close()


def get_dist_mat(values, test=spearman_roc, threads=8):
    """ Given a list of data, calculate the distances between them
        and return a NxN redundant array of these distances. """
    n = len(values)
    vals = []
    with Timer() as t:

        widgets = data.getWidgets('Evaluating distances: ')

        log("Creating {0}x{0} matrix, test function: ({1})"
            " {2} threads""".format(n, test.__name__, threads))

        # Generating matrix will be O(exp(n)) time
        c = list(combinations(values, 2))
        numcalc = len(c)

        pbar = pb.ProgressBar(widgets=widgets, maxval=numcalc)
        pbar.start()

        i = 0
        # Parallel code
        with concurrent.futures.ThreadPoolExecutor(8) as executor:
            batch = c[i:i+100]
            while batch:
                for val in executor.map(test, batch, timeout=30):
                    vals.append(val)
                    i += 1
                    pbar.update(i)
                batch = c[i:i+100]
        pbar.finish()

        vals = np.array(vals)
        mat = np.identity(n)
        """This step is key, basically we assign the upper triangle
          indices of a matrix size N, i.e.
          1    X    X    X
          0    1    X    X
          0    0    1    X
          0    0    0    1
          then we use the transpose of the matrix to copy the
          upper triangle into the lower triangle (making a
          symmetric matrix) """

        # Assign upper triangle
        try:
            mat[np.triu_indices(n, k=1)] = vals
        except ValueError as e:
            print("Error: {}".format(e))
            print("Couldn't broadcast array to triangle upper indices?")
            print("vals: {0}".format(vals))
            return

        # Make the matrix symmetric
        mat = (mat + mat.T) / 2
        np.fill_diagonal(mat, 0.0)

        # Because these tests give correlations not distances,
        # we must modify the values to give a distance equivalent
        if test is spearman_roc or test is kendall_tau:
            """ np.round() is used here because of floating point rounding
                (getting 1.0 - 1.0 != 0.0). Must perform this step to convert
                correlation data to distance """
            mat = 1.0 - np.round(mat, decimals=5)
            np.fill_diagonal(mat, 0.0)

        # Error checking for matrix to see if it is symmetric
        symmetry = np.allclose(mat.transpose(1, 0), mat)
        if not symmetry:
            logger.error('Matrix not symmetric...')

    log("Matrix took {0:.2}s to create: "
        "{1} pairwise calculations".format(t.elapsed(), numcalc))

    return mat


def cluster(mat, names, dump=None,
            dendrogram=None, method=None,
            distance=None):
    """ Takes an NxN array of distances and an array of names with
      the same indices, performs cluster analysis and shows a dendrogram"""

    clusters = []

    try:
        distArray = scipy.spatial.distance.squareform(mat)
    except ValueError as e:
        print(e)
        return clusters

    # This is the actual clustering using fastcluster
    Z = fc.linkage(distArray, method=method, metric=distance)
    clusters = scipy.cluster.hierarchy.fcluster(Z, 3., criterion='maxclust')

    if dendrogram:
        write_dendrogram_file(dendrogram, Z,
                              names,
                              distance=distance)
    return clusters


def key_from_indices(i, j, order=False):
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
        key = key_from_indices(i, j, order=order)
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


def bin_data_log(x, y, bins=10):
    return bin_data(np.log(x), np.log(y), bins)
