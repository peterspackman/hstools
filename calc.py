#!/usr/bin/python
import re
import sys
import numpy as np
import cv2
import matplotlib as mpl
from matplotlib import pyplot as plt
import scipy.stats as stats
import hist  as h
import scipy.spatial.distance
import fastcluster as fc
from scipy.cluster.hierarchy import dendrogram



def spearman_roc(hist1, hist2):
  """ Calculate the Spearman rank-order correlation coefficient from 2 histograms
  This may need to be modified, I'm uncertain whether or not 2 zeroes are ignored
  READ likely that they aren't, as such artificially high correlations are probable
  """ 
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
  if test==spearman_roc or test ==kendall_tau:
    mat = 1.0 - np.round(mat,decimals=5)
  return mat


def cluster(mat,names,tname):
  """ Takes an NxN array of distances and an array of names with
    the same indices, performs cluster analysis and shows a dendrogram
  """

  # Remove redundant distances from the NxN array (turning it into more of a triangle)
  distArray = scipy.spatial.distance.squareform(mat)
  
  # This is the clustering
  Z = fc.linkage(distArray,method='single',metric='euclidean')
  
  # Create a dendrogram
  R = dendrogram(Z,labels=names)

  # Plot stuff
  plt.xlabel('Compound Name')
  plt.ylabel('Dissimilarity')
  plt.suptitle('Clustering dendrogram of {0} compounds using {1}'.format(len(names),tname))
  plt.show()


