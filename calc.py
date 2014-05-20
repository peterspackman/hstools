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
def spearman_roc(hist1, hist2):
  x = hist1.flatten()
  y = hist2.flatten()
  return stats.spearmanr(x,y)

def kendall_tau(hist1, hist2):
  x = hist1.flatten()
  y = hist2.flatten()

  return stats.kendalltau(x,y)

#stats.kstest



