"""
lebedev.py

Contains wrapper methods to extract saved grids
from the h5file

"""
from os.path import dirname, abspath, join
import sbf
import numpy as np

DIR = dirname(abspath(__file__))
LEBEDEV_GRID_FILE = join(DIR, 'lebedev.sbf')
_GRIDS = sbf.read_file(LEBEDEV_GRID_FILE)

AVAILABLE_GRIDS = [i for i in range(3, 32, 2)] + [i for i in range(35, 132, 6)]
MAX_DEGREE = max(AVAILABLE_GRIDS)


def lebedev_grid(degree=21):
    """
    Returns the *angular* lebedev grid capable of exactly
    integrating a polynomial with given degree on the sphere.

    Grids are of shape(num_points, 3), with column 0
    being between (0, 2 Pi), column 1 between (0, Pi)
    and column 2 representing the weight for this grid point.

    >>> lebedev_grid(3)
    array([[3.14159265, 1.57079633, 0.16666667],
           [6.28318531, 1.57079633, 0.16666667],
           [4.71238898, 1.57079633, 0.16666667],
           [1.57079633, 1.57079633, 0.16666667],
           [4.71238898, 0.        , 0.16666667],
           [4.71238898, 3.14159265, 0.16666667]])
    """
    if degree > MAX_DEGREE:
        # error
        raise ValueError("maximum degree is {}".format(MAX_DEGREE))
    else:
        degree = next(x for x in AVAILABLE_GRIDS if x >= degree)
    rule = _GRIDS[str(degree)].data
    return rule


def integrate_lambda(grid, func):
    """
    Numerically integrate func f(theta, phi) on the quadrature given
    by grid.

    note that f must be applicable to numpy arrays using np.vectorize
    or similar

    Find the surface area of the unit sphere
    >>> import scipy.special
    >>> grid = lebedev_grid(31)
    >>> integrate_lambda(grid, lambda theta, phi: 1.0)
    12.566370614359423

    >>> l11 = lambda theta, phi: scipy.special.sph_harm(1, 1, theta, phi)
    >>> l11vec = np.vectorize(l11)
    >>> res = integrate_lambda(grid, l11vec)
    >>> all(np.isclose(res, [0,0]))
    True
    """

    return np.sum(func(grid[:, 0], grid[:, 1]) * grid[:, 2]) * 4 * np.pi


def integrate_values(grid, values):
    """
    Integrate the implicit function given by values associated with the grid.

    Find the surface of a unit sphere
    >>> grid = lebedev_grid(31)
    >>> integrate_values(grid, np.ones(len(grid)))
    12.566370614359425

    """
    return np.sum(values[:] * grid[:, 2] * 4 * np.pi)
