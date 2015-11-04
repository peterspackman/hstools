"""
    A module to store any global information
    such as the Van Der Waal's radii for various
    elements etc.
"""
import click
from timeit import default_timer

fg_colors = {'info':'white', 'error':'red', 'warning':'yellow'}

class Timer(object):
    """ A context manager timer class, to measure wall clock time"""
    def __init__(self):
        self.timer = default_timer
        self.start = 0
        self.end = 0
        self.elapsed_s = 0

    def __enter__(self):
        self.start = self.timer()
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.end = self.timer()
        self.elapsed_s = self.elapsed()
        self.elapsed_s = self.elapsed_s

    def elapsed(self):
        return self.timer() - self.start

    def __str__(self):
        if self.elapsed_s < 1.0:
            return '{:5.2f} ms'.format(self.elapsed_s * 1000.0)
        elif self.elapsed_s < 60.0:
            return '{:5.3f} s'.format(self.elapsed_s)
        else:
            secs = self.elapsed_s % 60.0
            mins = (self.elapsed_s - secs) // 60
            return '{} mins {:6.3f} s'.format(mins, secs)


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


def log(s, cat='info'):
    click.echo(click.style(s, fg=fg_colors[cat]))
