#!/usr/bin/python
# Core imports
from itertools import combinations
import json
import multiprocessing
import time

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
from data import log
import data
import pack.cio as cio


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
    """ Calculate a naive distance between two histograms"""
    H1, H2 = x
    hist1, _, _ = H1
    hist2, _, _ = H2
    x = hist1
    y = hist2
    dmat = x - y
    d = np.sum(dmat)
    return abs(d)


def dvalue(x):
    # unpack the tuple
    i1, i2 = x
    d = np.power(np.sum(np.power(i2 - i1, 2)), 0.5)
    return d


def get_dist_mat(histograms, test=spearman_roc, procs=4):
    """ Given a list of histograms, calculate the distances between them
        and return a NxN redundant array of these distances. """
    n = len(histograms)
    vals = []
    start_time = time.time()
    widgets = data.widgets

    output = "Creating {0}x{0} matrix, test={1}, using {2} processes"
    log(output.format(n, test.__name__, procs))

    # Here is why our asymptote is exponential! :(
    c = list(combinations(histograms, 2))
    numcalc = len(c)

    pbar = pb.ProgressBar(widgets=widgets, maxval=numcalc)
    pbar.start()
    # Parallel code
    p = multiprocessing.Pool(procs)

    r = p.map_async(test, c, callback=vals.extend)
    p.close()
    while True:
        if r.ready():
            break
        pbar.update()
        time.sleep(0.2)
    p.join()
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
    mat[np.triu_indices(n, k=1)] = vals
    # Make the matrix symmetric
    mat = (mat + mat.T) / 2

    # Because these tests give correlations not distances,
    # we must modify the values to give a distance equivalent

    if test is spearman_roc or test is kendall_tau:
        """ np.round() is used here because of floating point rounding
            (getting 1.0 - 1.0 != 0.0). Must perform this step to convert
            correlation data to distance """
        mat = 1.0 - np.round(mat, decimals=5)
    else:
        np.fill_diagonal(mat, 0.0)
    t = time.time() - start_time
    output = 'Matrix took {0:.2}s to create, performing {1} calculations'
    log(output.format(t, numcalc))
    return mat


def cluster(mat, names, tname, dump=None,
            dendrogram=None, method=None,
            distance=None):
    """ Takes an NxN array of distances and an array of names with
      the same indices, performs cluster analysis and shows a dendrogram"""

    distArray = scipy.spatial.distance.squareform(mat)
    start_time = time.time()

    # This is the actual clustering using fastcluster
    Z = fc.linkage(distArray, method=method, metric=distance)
    outstring = 'Clustering {0} histograms'.format(len(names))
    outstring += ' took {0:.3}s'.format(time.time() - start_time)
    log(outstring)
    threshold = distance*max(Z[:, 2])
    if dendrogram:
        # Create a dendrogram
        dend(Z, labels=names, color_threshold=threshold)
        # Plot stuff
        plt.xlabel('Compound Name')
        plt.ylabel('Dissimilarity')
        dpi = 200
        plt.suptitle("""Clustering dendrogram of {0}
                    compounds using {1}""".format(len(names), tname))
        if len(names) > 100:
            fig = plt.gcf()
            fig.set_size_inches(10.5, min(len(names)*0.1, 32768/dpi))
        log('Saving dendrogram')
        plt.savefig(dendrogram, dpi=dpi)
        plt.close()

    if dump:
        log('Dumping tree structure in {0}'.format(dump))
        T = scipy.cluster.hierarchy.to_tree(Z, rd=False)
        d = dict(children=[], name="Root1")
        add_node(T, d)
        label_tree(d["children"][0], names)
        json.dump(d, open(dump, 'w'), sort_keys=True, indent=4)
        log('printing clusters')
        # HARDCODED NUMBER OF CLUSTERS
        clusters = scipy.cluster.hierarchy.fcluster(Z, 10,
                                                    criterion='maxclust')
        nclusters = clusters.size
        num = max(clusters)
        c = []
        for i in range(1, num):
            c.append([names[x] for x in range(nclusters) if clusters[x] == i])
        json.dump(c, open('clusters.txt', 'w'), indent=4)


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
    """ Helper function to label the tree """
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

    n["name"] = "-".join(sorted(map(str, leafNames)))
    if len(n["name"]) > 16:
        n["name"] = n["name"][:16] + '...'
    # Labeling convention: "-" separates leaf names

    return leafNames


def area_tri(a, b, c):
    """ Calculate the area of a triangle given by its 3 vertices
    using the cross product formula |AxB|/2"""
    v1 = a - b
    v2 = c - b
    # Because these functions expect 6 doubles or 3 doubles as arguments
    # we have to unpack the values!
    x, y, z = cio.cross3D(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2])
    return cio.normal3D(x, y, z) / 2


def get_contrib_percentage(vertices, indices, internal,
                           external, distances,
                           dp=8, restrict=True,
                           order=False):
    """ Given a the triangles that make up a hirshfeld surface,
    and lists of the closest internal and external atoms along
    with their respective distances from the surface,
    calculate the makeup of the hirshfeld surface in terms of
    which element->element interactions are responsible for that
    area """
    contrib = {}
    contrib_p = {}
    if restrict:  # Check if we can restrict
        unique = np.unique(np.append(internal, external))
        for sym in unique:
            if sym not in data.vdw_radii:
                restrict = False

    for i in range(internal.size):
        # Key in the form "internal -> external" e.g. "F -> H"
        chsymi = internal[i]
        chsyme = external[i]
        if(not order):
            chsymi, chsyme = sorted((chsymi, chsyme))
        key = "{0} -> {1}".format(chsymi, chsyme)
        if restrict:
            avg_d = np.mean(distances[indices[i]])
            threshold = data.vdw_radii[chsymi] + data.vdw_radii[chsyme]
        if not restrict or avg_d < threshold:
            tri = [vertices[n] for n in indices[i]]
            area = area_tri(tri[0], tri[1], tri[2])
            if key in contrib:
                contrib[key] += area
            else:
                contrib[key] = area

    for x in contrib:
        p = np.round(contrib[x] / sum(contrib.values()), decimals=8)
        contrib_p[x] = p

    return contrib, contrib_p
