#!/usr/bin/python
# Core imports
import re
import sys
import time
import json
from itertools import combinations
import multiprocessing

# Library imports
import matplotlib as mpl
from matplotlib import pyplot as plt
import scipy.spatial.distance
import fastcluster as fc
import numpy as np
import scipy.cluster.hierarchy
from scipy.cluster.hierarchy import dendrogram
import scipy.stats as stats
import progressbar as pb

# Local imports
import hist as h
import pack.cio as cio
from data import vdw_radii


def spearman_roc(x):
    """ Calculate the Spearman rank-order correlation coefficient from
    2 histograms This may need to be modified, I'm uncertain whether or
    not 2 zeroes are ignored READ likely that they aren't, as such
    artificially high correlations are probable"""
    H1, H2 = x
    hist1, _, _ = H1
    hist2, _, _ = H2
    x = hist1.flatten()
    y = hist2.flatten()

    r, p = stats.spearmanr(x, y)
    return r


def kendall_tau(x):
    """ Calculate Kendall's Tau from the given histograms"""
    H1, H2 = x
    hist1, _, _ = H1
    hist2, _, _ = H2
    x = hist1.flatten()
    y = hist2.flatten()
    r, p = stats.kendalltau(x, y)
    return r


def hdistance(x):
    H1, H2 = x
    hist1, _, _ = H1
    hist2, _, _ = H2
    x = hist1
    y = hist2
    dmat = x - y
    d = np.sum(dmat)
    return abs(d)


def get_dist_mat(histograms, test=spearman_roc, processes=4):
    """
        Given a list of histograms, calculate the distances between them
        and return a NxN redundant array of these distances

        It should be noted that there is a potential inefficiency
        here as x-> should be the same as y->x, so we could cut cpu time in
        half roughly but cutting out the inefficiency."""
    n = len(histograms)

    output = "Creating {0}x{0} matrix, test={1}, using {2} processes"
    print output.format(n, test.__name__, processes)
    start_time = time.time()
    widgets = ['Reading Files: ', pb.Percentage(), ' ',
               pb.Bar(marker='=', left='[', right=']'),
               ' ', pb.ETA()]
    c = list(combinations(histograms, 2))
    numcalc = len(c)
    pbar = pb.ProgressBar(widgets=widgets, maxval=numcalc)
    p = multiprocessing.Pool(processes)
    vals = []

    r = p.map_async(test, c, callback=vals.extend)
    p.close()
    done = 0
    while True:
        if r.ready():
            break
        if numcalc - r._number_left > done:
            pbar.update(done)
            done = numcalc - r._number_left
        time.sleep(0.2)
    p.join()
    pbar.finish()

    # vals = map(test, combinations(histograms, 2))
    vals = np.array(vals)
    mat = np.identity(n)
    """
      This step is key, basically we select upper triangle
      indices of size N, i.e.
      0    1    1    1
      0    0    1    1
      0    0    0    1
      0    0    0    0
    """

    mat[np.triu_indices(n, k=1)] = vals
    # Make the matrix symmetric
    mat = (mat + mat.T) / 2

    # Because these tests give correlations not distances,
    # we must modify the values to give a distance equivalent

    if test is spearman_roc or test is kendall_tau:
        # np.round() is used here because of floating point rounding
        # (getting 1.0 - 1.0 != 0.0)
        mat = 1.0 - np.round(mat, decimals=5)
    t = round(time.time() - start_time, 3)
    print 'matrix took {0:.2} seconds to create'.format(t)
    return mat


def cluster(mat, names, tname, dump=None):
    """ Takes an NxN array of distances and an array of names with
      the same indices, performs cluster analysis and shows a dendrogram
    """

    distArray = scipy.spatial.distance.squareform(mat)
    start_time = time.time()

    # This is the actual clustering using fastcluster
    Z = fc.linkage(distArray, method='single', metric='euclidean')
    outstring = 'Clustering {0} histograms'.format(len(names))
    outstring += ' took {0:.3}s'.format(time.time() - start_time)
    print outstring
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
        print 'Dumping tree structure in {0}'.format(dump)
        T = scipy.cluster.hierarchy.to_tree(Z, rd=False)
        d = dict(children=[], name="Root1")
        add_node(T, d)
        label_tree(d["children"][0], names)
        json.dump(d, open(dump, 'w'), sort_keys=True, indent=4)


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
    v1 = a - b
    v2 = c - b
    # Because these funcitons expect 6 doubles or 3 doubles as arguments
    # we have to unpack the values!
    x, y, z = cio.cross3D(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2])
    return cio.normal3D(x, y, z) / 2


def get_contrib_percentage(vertices, indices, internal, external, dp=3):
    contrib = {}
    contrib_p = {}
    RESTRICT_DISTANCE = False

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
