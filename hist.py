#!/usr/bin/python
# This module contains methods for the creation and manipulation of histograms
import re
import sys
import numpy as np
import scipy.sparse

def bin_data(x, y, bins = 10):
  nx = ny = bins
  xmin, xmax = min(x), max(x)
  ymin, ymax = min(y), max(y)
  dx = (xmax -xmin) / (nx - 1.0)
  dy = (ymax - ymin) / (ny -1.0)

  weights = np.ones(len(x))

  # this is a slightly modified version of np.digitize()

  xyi = np.vstack((x,y)).T
  xyi -= [xmin,ymin]
  xyi /= [dx, dy]
  xyi = np.floor(xyi, xyi).T

  grid = scipy.sparse.coo_matrix((weights, xyi), shape=(nx, ny)).toarray()

  return grid, np.linspace(xmin, xmax, nx), np.linspace(ymin,ymax,ny)

def bin_data_log(x, y, bins = 10):
  return bin_data(np.log(x), np.log(y), bins)
  


  


