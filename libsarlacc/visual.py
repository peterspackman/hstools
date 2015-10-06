#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans, vq

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']


def scatter(x, y):
    """
    Produce a 2D scatter plot of x and y
    """
    data = np.array(zip(x, y))
    n = 7
    centroids, _ = kmeans(data, n)
    idx, _ = vq(data, centroids)

    for i in range(0, n):
        plt.scatter(data[idx == i, 0], data[idx == i, 1], c=colors[i])

    plt.plot(centroids[:, 0], centroids[:, 1], 'sg', markersize=8)
    plt.show()
