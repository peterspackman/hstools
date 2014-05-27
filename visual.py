#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans,vq
colors = ['b','g','r','c','m','y','k']

def scatter(x,y,km=True):
  """
  Produce a 2D scatter plot of x and y
  """
  data = np.array(zip(x,y))
  print data
  n = 7
  centroids, _ = kmeans(data,n)
  idx, _ = vq(data,centroids)


  for i in range(0,n):
    plt.scatter(data[idx==i,0],data[idx==i,1],c=colors[i])
  plt.plot(centroids[:,0],centroids[:,1],'sg',markersize=8)
  plt.show()


class DripHist:
  "A class to construct a drip histogram"

  def __init__(self,radius=0.05):
    UserDict.__init__(self)
    self.drops = []
 
  def getDrops(self):
    return self.drops


  def add_point(p):
    pass
