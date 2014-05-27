#!/usr/bin/python
import re
import sys
import numpy as np
import cv2
import matplotlib as mpl
from matplotlib import pyplot as plt
import scipy.stats as stats
import hist  as h


# Calculate the Spearman rank-order correlation coefficient from 2 histograms
# This may need to be modified, I'm uncertain whether or not 2 zeroes are ignored
# READ likely that they aren't, as such artificially high correlations are probable

def spearman_roc(hist1, hist2):
  x = hist1.flatten()
  y = hist2.flatten()
  return stats.spearmanr(x,y)

def kendall_tau(hist1, hist2):
  x = hist1.flatten()
  y = hist2.flatten()

  return stats.kendalltau(x,y)

#stats.kstest


def get_correl_mat(histograms,test=spearman_roc):
  n = len(histograms)

  print "Creating a {0}x{0} matrix using coefficients from {1}".format(n,test.__name__)

  mat = np.zeros( (n,n) )
  
  for i, (H1,x1,y1) in enumerate(histograms):
    for j, (H2,x2,y2) in enumerate(histograms):
      r, p = test(np.array(H1),np.array(H2))
      mat[i][j] = r
  return mat



