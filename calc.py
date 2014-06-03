#!/usr/bin/python
# Core imports
import re
import sys
import time
# Library imports
import matplotlib as mpl
from matplotlib import pyplot as plt
import scipy.spatial.distance
import fastcluster as fc
import numpy as np
from scipy.cluster.hierarchy import dendrogram
import scipy.stats as stats
# Local imports
import hist as h


def spearman_roc(H1, H2):
    """ Calculate the Spearman rank-order correlation coefficient from
    2 histograms This may need to be modified, I'm uncertain whether or
    not 2 zeroes are ignored READ likely that they aren't, as such
    artificially high correlations are probable"""
    hist1, _, _ = H1
    hist2, _, _ = H2
    x = hist1.flatten()
    y = hist2.flatten()

    r, p = stats.spearmanr(x, y)
    return r


def kendall_tau(H1, H2):
    """ Calculate Kendall's Tau from the given histograms"""
    hist1, _, _ = H1
    hist2, _, _ = H2
    x = hist1.flatten()
    y = hist2.flatten()
    r, p = stats.kendalltau(x, y)
    return r


def get_correl_mat(histograms, test=spearman_roc):
    """ SHOULD BE RENAMED TO get_dist_mat
        Given a list of histograms, calculate the distances between them
        and return a NxN redundant array of these distances

        It should be noted that there is a potential inefficiency
        here as x-> should be the same as y->x, so we could cut cpu time in
        half roughly but cutting out the inefficiency."""
    n = len(histograms)

    output = "Creating a {0}x{0} matrix using coefficients from {1}"
    print output.format(n, test.__name__)
    start_time = time.time()

    # this will be O(n^2) but I don't see a way around that
    mat = [[test(H1, H2) for H1 in histograms] for H2 in histograms]

    # Ensure it is a np.array (though it should be already)
    mat = np.array(mat)

    # Because these tests give correlations not distances,
    # we must modify the values to give a distance equivalent

    if test is spearman_roc or test is kendall_tau:
        # np.round() is used here because of floating point rounding
        # (getting 1.0 - 1.0 != 0.0)
        mat = 1.0 - np.round(mat, decimals=5)

    print 'matrix took {0} seconds to create'.format(time.time() - start_time)
    return mat


def cluster(mat, names, tname):
    """ Takes an NxN array of distances and an array of names with
      the same indices, performs cluster analysis and shows a dendrogram
    """

    distArray = scipy.spatial.distance.squareform(mat)
    start_time = time.time()

    # This is the actual clustering using fastcluster
    print 'Clustering {0} data points...'.format(len(names))
    Z = fc.linkage(distArray, method='single', metric='euclidean')
    print 'took {0} seconds'.format(time.time() - start_time)
    # Create a dendrogram
    R = dendrogram(Z, labels=names)

    # Plot stuff
    plt.xlabel('Compound Name')
    plt.ylabel('Dissimilarity')
    plt.suptitle("""Clustering dendrogram of {0}
                 compounds using {1}""".format(len(names), tname))
    plt.savefig('dendrogram.png')
    plt.close()


def get_contrib_percentage(internal, external, dp=3):
    contrib = {}
    contrib_p = {}

    for i in range(internal.size):
        # Key in the form "internal -> external" e.g. "F -> H"
        key = "{0} -> {1}".format(internal[i], external[i])
        if key in contrib:
            contrib[key] += 1
        else:
            contrib[key] = 1

    for x in contrib:
        p = np.round(contrib[x] * 100.0 / sum(contrib.values()), decimals=dp)
        contrib_p[x] = p

    return contrib, contrib_p
