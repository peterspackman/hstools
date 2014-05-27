#!/usr/bin/python
import sys
import scipy
import hist
import numpy as np
import fileio as fio
import cv2 as cv

def num_hotspots(h,n,threshold = 0.05):
  h = h.flatten()
  vals = np.where(h >n*threshold,1,0)
  print vals
  return sum(vals)

def compute_histogram(src,h_bins=30,s_bins=32,scale=10):
  #create images
  src = cv.cvtColor(src,cv.CV_BGR2HSV)

  hrange = (0,180)
  srange = (0,256)
  ranges = [hrange,srange]

  channels = [0,1]

  
  #compute histogram
  hist = cv.calcHist(src,channels,cv.Mat(),
  return hist

def main():
  src_base = cv.imread(sys.argv[1])
  src_test1 = cv.imread(sys.argv[2])
  src_test2 =cv.imread(sys.argv[3])

  hist_base = compute_histogram(src_base)
  hist1 = compute_histogram(src_test1)
  hist2 = compute_histogram(src_test2)

  sc = cv.CompareHist(hist1,hist2,cv.CV_COMP_BHATTACHARYYA)

  print sc

  planes = [hplane,splane]

  


if __name__ == '__main__':
 main()

