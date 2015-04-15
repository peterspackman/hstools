"""
    A module to store any global information
    such as the Van Der Waal's radii for various
    elements etc.
"""
import logging
import progressbar as pb
from timeit import default_timer
from periodictable import elements

FORMAT = "%(levelname)s%(message)s"
logging.basicConfig(level=logging.INFO, format=FORMAT)
logging.addLevelName(logging.ERROR, 'error: ')
logging.addLevelName(logging.WARNING, 'warning: ')
logging.addLevelName(logging.INFO, '')
logging.addLevelName(logging.DEBUG, 'debug: ')
logging.addLevelName(logging.CRITICAL, 'CRITICAL: ')

logger = logging.getLogger("sarlacc")


class Timer(object):
    """ A context manager timer class, to measure wall clock time"""
    def __init__(self):
        self.timer = default_timer

    def __enter__(self):
        self.start = self.timer()
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.end = self.timer()
        self.elapsed_s = self.elapsed()
        self.elapsed_s = self.elapsed_s * 1000.0

    def elapsed(self):
        return self.timer() - self.start


def log_traceback(e):
    logger.exception(e)


def logClosestPair(mat, names):
    import numpy as np
    np.fill_diagonal(mat, np.inf)
    x = np.nanargmin(mat)
    ind = (x//len(names), x % len(names))
    a = str(names[ind[0]])
    b = str(names[ind[1]])
    log('Closest pair: {0}, d= {1:.5f}'.format((a, b), mat[ind]))
    np.fill_diagonal(mat, 0.0)


def logFarthestPair(mat, names):
    import numpy as np
    x = np.nanargmax(mat)
    ind = (x//len(names), x % len(names))
    a = str(names[ind[0]])
    b = str(names[ind[1]])
    log('Farthest pair: {0}, d= {1:.5f}'.format((a, b), mat[ind]))


def getWidgets(msg, color='white'):
    return [msg, pb.Percentage(), ' ',
            pb.Bar(marker=chr(0x2500), left='',
            right=''), ' ', pb.ETA(), ' ']


def log(s):
    logger.info(s)


def log_error(s):
    logger.error(s)
