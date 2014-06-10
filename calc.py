#!/usr/bin/python
# Core imports
import re
import sys
import time
import json
# Library imports
import matplotlib as mpl
from matplotlib import pyplot as plt
import scipy.spatial.distance
import fastcluster as fc
import numpy as np
import scipy.cluster.hierarchy
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


def hdistance(H1, H2):
    hist1, _, _ = H1
    hist2, _, _ = H2
    x = hist1
    y = hist2
    dmat = x - y
    d = np.sum(dmat)
    return abs(d)


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


def cluster(mat, names, tname, dump=False):
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
    plt.savefig('dendrogram.png', dpi=800)
    plt.close()
    if dump:
        d3 = "d3-dendrogram.json"
        print 'Dumping tree structure in {0}'.format(d3)
        T = scipy.cluster.hierarchy.to_tree(Z, rd=False)
        d = dict(children=[], name="Root1")
        add_node(T, d)
        label_tree(d["children"][0], names)
        json.dump(d, open(d3, 'w'), sort_keys=True, indent=4)


def add_node(node, parent):
    """ A Helper method for outputting the dendrogram
    linkage for visualisation in d3.js"""
    newNode = dict(node_id=node.id, children=[])
    parent["children"].append(newNode)
    # Recursively add the current node's children
    if node.left:
        add_node(node.left, newNode)
    if node.right:
        add_node(node.right, newNode)


def label_tree(n, names):
    id2name = dict(zip(range(len(names)), names))
    # If it's a leaf node we have the name
    if len(n["children"]) == 0:
        leafNames = [id2name[n["node_id"]]]
    # Otherwise flatten all the leaves in the subtree
    else:
        leafNames = reduce(lambda ls, c: ls + label_tree(c, names),
                           n["children"], [])
    # Delete the node id as it is no longer needed
    del n["node_id"]

    n["name"] = name = "-".join(sorted(map(str, leafNames)))
    if len(n["name"]) > 16:
        n["name"] = n["name"][:16] + '...'
    # Labeling convention: "-" separates leaf names

    return leafNames


# Calculate the area of a triangle with given by points a,b,c
def area_tri(a, b, c):
    return np.linalg.norm(np.cross(a - b, c - b)) / 2


def get_contrib_percentage(vertices, indices, internal, external, dp=3):
    contrib = {}
    contrib_p = {}

    for i in range(internal.size):
        # Key in the form "internal -> external" e.g. "F -> H"
        key = "{0} -> {1}".format(internal[i], external[i])
        tri = [vertices[n] for n in indices[i]]
        area = area_tri(tri[0], tri[1], tri[2])
        if key in contrib:
            contrib[key] += area
        else:
            contrib[key] = area

    for x in contrib:
        p = np.round(contrib[x] * 100.0 / sum(contrib.values()), decimals=dp)
        contrib_p[x] = p

    return contrib, contrib_p
