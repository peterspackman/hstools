"""
    A module to store any global information
    or configuration
"""
from timeit import default_timer

PREFIX = {'info': '', 'error': 'ERR: ', 'warning': 'warning: '}


class Timer(object):
    """ A context manager timer class, to measure wall clock time"""
    def __init__(self):
        self.timer = default_timer
        self.begin = 0
        self.end = 0
        self.elapsed_s = 0

    def __enter__(self):
        self.begin = self.timer()
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.end = self.timer()
        self.elapsed_s = self.elapsed()
        self.elapsed_s = self.elapsed_s

    def start(self):
        """ Start timing! """
        self.begin = self.timer()

    def elapsed(self):
        """ Return the time difference between now and when
        the object was entered"""
        return self.timer() - self.begin

    def __str__(self):
        if self.elapsed_s < 1.0:
            return '{:5.2f} ms'.format(self.elapsed_s * 1000.0)
        elif self.elapsed_s < 60.0:
            return '{:5.3f} s'.format(self.elapsed_s)
        else:
            secs = self.elapsed_s % 60.0
            mins = (self.elapsed_s - secs) // 60
            return '{} mins {:6.3f} s'.format(mins, secs)

def log(string, cat='info'):
    """Write s to stdout, with prefix based on category"""
    print(PREFIX[cat] + string)
