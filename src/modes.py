# Core imports
import sys
import time
# Library imports
import numpy as np
# Local imports
from .data import log, log_traceback


def logClosestPair(mat, names):
    np.fill_diagonal(mat, np.inf)
    x = np.nanargmin(mat)
    ind = (x//len(names), x % len(names))
    a = str(names[ind[0]])
    b = str(names[ind[1]])
    log('Closest pair: {0}, d= {1:.5f}'.format((a, b), mat[ind]))
    np.fill_diagonal(mat, 0.0)


def logFarthestPair(mat, names):
    x = np.nanargmax(mat)
    ind = (x//len(names), x % len(names))
    a = str(names[ind[0]])
    b = str(names[ind[1]])
    log('Farthest pair: {0}, d= {1:.5f}'.format((a, b), mat[ind]))
