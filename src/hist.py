#!/usr/bin/python
""" This module contains methods for the creation and
    manipulation of histograms from given data"""
# Core imports
import re
import sys
# Library imports
import numpy as np
import scipy.sparse


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
