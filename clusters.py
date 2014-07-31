from __future__ import division
import json
import sys


def check_clusters(file1, file2):
    clusters1 = json.load(open(file1, 'r'))
    clusters2 = json.load(open(file2, 'r'))
    for cluster1 in clusters1:
        for cluster2 in clusters2:
            n = max(len(cluster1), len(cluster2))
            intersect = len(set(cluster1).intersection(cluster2))
            overlap = intersect / n
            if overlap > 0.0:
                print '{0:.0%} overlap between clusters'.format(intersect / n)


if __name__ == '__main__':
    check_clusters(sys.argv[1], sys.argv[2])
