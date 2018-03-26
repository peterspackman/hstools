"""Contains logic associated with a spherical harmonic transform"""
import logging
import numpy as np
from scipy.special import sph_harm
from .lebedev import lebedev_grid, integrate_values

LOG = logging.getLogger(__name__)

class SHT:
    """Encapsulate logic of spherical harmonic transform implementations"""
    _shtns = None
    _l_max = 2
    def __init__(self, l_max, force_python=False):
        """initialize a spherical harmonic transform object"""
        self._l_max = l_max
        if not force_python:
            try:
                import shtns
                self._shtns = shtns.sht(l_max, l_max)
            except ImportError:
                LOG.debug('Could not import shtns')

    @property
    def grid(self):
        """Return the set of angular grid points for this sht"""
        if self._shtns:
            _, nphi = self._shtns.set_grid()
            phi, theta = np.meshgrid(np.arccos(self._shtns.cos_theta),
                                     np.arange(nphi)*(2*np.pi/nphi))
            return np.vstack((theta.flatten(), phi.flatten())).transpose()
        return lebedev_grid(degree=131)

    def analyse(self, values):
        """Perform a spherical harmonic transform given a grid and a set of values

        Arguments:
        values -- set of complex scalar function values associated with grid points
        """
        if self._shtns:
            desired_shape = self._shtns.spat_shape[::-1]
            return self._shtns.analys_cplx(
                np.array(values, dtype=np.complex128).reshape(desired_shape).transpose())
        coefficients = np.zeros((self.l_max + 1)*(self.l_max + 1),
                                dtype=np.complex128)
        lm = 0
        grid = self.grid
        for l in range(0, self.l_max + 1):
            for m in range(-l, l + 1):
                vals = np.conj(sph_harm(m, l, grid[:, 0], grid[:, 1]))
                vals *= values
                coefficients[lm] = integrate_values(grid, vals)
                lm += 1
        return coefficients


    def synthesis(self, coefficients):
        """Perform a spherical harmonic transform given a grid and a set of values

        Arguments:
        values -- set of complex scalar function values associated with grid points
        """
        if self._shtns:
            max_coeff = (self.l_max+1)**2
            return self._shtns.synth_cplx(coefficients[:max_coeff]).transpose().flatten()

        # Fall back to the (much slower) python implementation
        grid = self.grid
        values = np.zeros(len(grid), dtype=np.complex)
        l_m = 0
        for l in range(0, self.l_max + 1):
            for m in range(-l, l + 1):
                ylm = sph_harm(m, l, grid[:, 0], grid[:, 1])
                values[:] += coefficients[l_m] * ylm
                l_m += 1
        return values


    @property
    def l_max(self):
        """Maximum angular momenta to evaluate up to"""
        return self._l_max
