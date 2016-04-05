# Core imports
from collections import defaultdict, namedtuple
import concurrent.futures
import glob
import os
import sys
from functools import reduce

# Library imports
import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Local imports
from .calc import bin_data, get_contrib
from .config import (
        Timer,
        log,
)

from pathlib import Path


class lazy_property(object):
    """
    Helper class to be used for lazy evaluation of an object attribute.
    property should represent non-mutable data, as it replaces itself.
    """
    def __init__(self, fget):
        self.fget = fget
        self.func_name = fget.__name__

    def __get__(self, obj, cls):
        if obj is None:
            return None
        value = self.fget(obj)
        setattr(obj, self.func_name, value)
        return value

HarmonicsData = namedtuple('HarmonicsData',
                           'radius coefficients invariants name')
_SurfaceDataTuple = namedtuple('SurfaceData',
                               'contributions formula name')
_FingerprintDataTuple = namedtuple('FingerprintData',
                                   'd_e d_i name')


class FingerprintData(_FingerprintDataTuple):
    """
    Light wrapper around namedtuple for
    the one-time calculation of the histogram
    from given $d_e$ and $d_i$ values on
    the Hirshfeld surface.
    """
    @lazy_property
    def histogram(self):
        h = bin_data(self.d_e, self.d_i, bins=6)
        return h


class SurfaceData(_SurfaceDataTuple):
    """
    Light wrapper around namedtuple for
    the one-time calculation of the surface area contributions
    from given $d_e$ and $d_i$ values on
    the Hirshfeld surface.
    """

    @lazy_property
    def contributions_dict(self):
        _, c = get_contrib(self.contributions)
        return c


class DataFileReader:
    """ Read hdf5 data files, given a reader object"""
    def __init__(self, attributes, obj):
        self.attributes = attributes
        self.object = obj

    def read(self, path):
        outputs = []
        try:
            with h5py.File(str(path), 'r') as f:
                for group in f:
                    output = {}
                    for k, v in self.attributes.items():
                        if v not in f[group]:
                            log("Couldn't find dataset: {} in group {}".format(v, g),
                                cat='error')
                        output[k] = np.array(f[group][v])
                    output['name'] = path.stem + ('-' + group)
                    outputs.append(output)
            return [self.object(**output) for output in outputs]

        except TypeError as e:
            log(str(e))
        except Exception as e:
            log("Ignoring '{}'; Caught ({})".format(path.name,
                                                    type(e).__name__),
                cat='warning')
            log(str(e))
            pass


def load_saved_data(fname):
    attributes = ['matrix', 'names']
    data = readh5file(fname, attributes)
    return data['matrix'], data['names']


# Takes a dictionary of dataset_name
def write_hdf5_file(fname, attributes):
    with h5py.File(fname, 'w') as f:
        for k in attributes.keys():
            data = attributes[k]
            dset = f.create_dataset(k,
                                    data.shape,
                                    data.dtype,
                                    compression='gzip')
            dset[...] = data[...]


def standard_figure(extent=[0.5, 2.5, 0.5, 2.5], figsize=(4, 4), dpi=300):
    plt.style.use('seaborn-white')
    f = plt.figure(figsize=figsize, dpi=dpi)
    cmap = mpl.cm.viridis
    ax = f.add_subplot(111, xlim=extent[0:2], ylim=extent[2:4])
    plt.xticks(np.arange(extent[0], extent[1], (extent[1] - extent[0])/5))
    plt.yticks(np.arange(extent[2], extent[3], (extent[3] - extent[2])/5))
    plt.grid(b=True, which='major', axis='both')

    ticklabelpad = mpl.rcParams['xtick.major.pad']

    # Label our graph axes in nice places!
    plt.annotate(r'$d_i$', fontsize=14, xy=(1, 0), xytext=(5, - ticklabelpad),
                 ha='left', va='top', xycoords='axes fraction',
                 textcoords='offset points')
    plt.annotate(r'$d_e$', fontsize=14, xy=(0, 1.02),
                 xytext=(5, - ticklabelpad), ha='right',
                 va='bottom', xycoords='axes fraction',
                 textcoords='offset points')
    return f, cmap, ax, extent


def kde_plotfile(x, y, fname='outhex.png'):
    import seaborn as sns
    f, cmap, ax, _ = standard_figure()

    sns.kdeplot(x, y, cmap=cmap, shade=True, ax=ax)
    ax.collections[0].set_alpha(0)

    f.savefig(fname, bbox_inches='tight')
    plt.clf()
    plt.close()


def hexbin_plotfile(x, y, fname='outhex.png', kind='log', nbins=100):
    extent = [min(0.5, np.min(x)), max(2.5, np.max(x)*1.2),
              min(0.5, np.min(y)), max(2.5, np.max(y)*1.2)]
    f, cmap, _, extent = standard_figure(extent=extent)

    plt.hexbin(x, y, mincnt=1, gridsize=nbins,
               cmap=cmap, bins=kind,
               extent=extent)

    f.savefig(fname, bbox_inches='tight')
    plt.clf()
    plt.close()


# FILE FUNCTIONS
def plotfile(x, y, fname='out.png', kind='linear', nbins=100):
    """ Construct a histogram, plot it, then write the image as a png
    to a given filename."""

    # Not sure why but we have linear and log as options
    if(kind == 'linear'):
        H, xedges, yedges = bin_data(x, y, bins=nbins)

    # NEED TO ROTATE AXES
    H = np.rot90(H)
    H = np.flipud(H)

    f, cmap, ax, _ = standard_figure()

    # this is the key step, the rest is formatting
    plt.pcolormesh(xedges, yedges, H,
                   cmap=cmap,
                   norm=mpl.colors.LogNorm())

    plt.grid(b=True, which='major', axis='both')

    # SAVE TO PNG WITH LITTLE TO NO BORDER
    f.savefig(fname, bbox_inches='tight')

    # CLOSE UP
    plt.clf()
    plt.close()


def write_sa_file(fname, surfaces):
    """ Write out the surface area contribution info to a given
    file."""
    with open(fname, 'w') as f:
        for x in surfaces:
            contributions = x.contributions
            line = '{0}, {1}'.format(x.name, x.formula)
            if not contributions:
                line = line + '--Nil--'
            else:
                for key in sorted(contributions,
                                  key=lambda key: contributions[key]):
                    line = line + ', '
                    line = line + '{} = {:.2%}'.format(key, contributions[key])
            f.write(line + '\n')


def write_mat_file(fname, mat, names, clusters):
    d = defaultdict()
    log("Writing clustering data to  '{}'".format(str(fname)))
    d['matrix'] = np.array(mat)
    d['names'] = np.array(names, dtype='<S64')
    d['clusters'] = np.array(clusters)
    write_hdf5_file(fname, d)


def batch_process(files, reader, procs=4):
    """
    Takes a list of files and a reader class, and
    constructs instances of the given class from data
    in the files using ProcessPoolExecutor
    """
    n = len(files)
    with Timer() as t:
        vals = []

        with concurrent.futures.ProcessPoolExecutor(procs) as executor:
            fs = [executor.submit(reader.read, f) for f in files]
            for i, f in enumerate(concurrent.futures.as_completed(fs)):
                vals.append(f.result())

        vals = [x for sublist in vals for x in sublist if sublist is not None]

        n_success = len(vals)
        if n_success < n:
            err = 'Skipped {0} files due to errors'.format(n - n_success)
            log(err, cat='error')

    log("Read {1} files in {0}".format(t, n))

    return sorted(vals, key=lambda x: x.name)
