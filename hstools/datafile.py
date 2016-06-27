"""
    Collection of classes and methods related to reading HS
    data from HDF5 files
"""
# Core imports
from collections import defaultdict, namedtuple
import concurrent.futures

# Library imports
import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import viridis

# Local imports
from .calc import bin_data, get_contrib
from .config import Timer, log


def lazy_property(func):
    """
    Helper method to be used for lazy evaluation of an object attribute.
    """
    attr_name = '_lazy_' + func.__name__
    @property
    def _lazy_property(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, func(self))
        return getattr(self, attr_name)
    return _lazy_property

HarmonicsData = namedtuple('HarmonicsData',
                           'radius coefficients invariants name')
_SurfaceDataTuple = namedtuple('SurfaceData',
                               'contributions name')
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
        """ Calculate the histogram from raw d_e/d_i values"""
        return bin_data(self.d_e, self.d_i, bins=6)


class SurfaceData(_SurfaceDataTuple):
    """
    Light wrapper around namedtuple for
    the one-time calculation of the surface area contributions
    from given $d_e$ and $d_i$ values on
    the Hirshfeld surface.
    """

    @lazy_property
    def contributions_dict(self):
        """Retrieve a dictionary of surface area contributions"""
        _, contribution_percentage = get_contrib(self.contributions)
        return contribution_percentage


class DataFileReader:
    """ Read hdf5 data files, given a reader object"""
    def __init__(self, attributes, obj):
        self.attributes = attributes
        self.object = obj

    def read(self, path):
        """Given a HDF5 file, reads relevant data from each group and returns a
            list of objects of kind self.obj"""
        outputs = []
        try:
            with h5py.File(str(path), 'r') as root:
                for group in root:
                    output = {}
                    for key, dataset_name in self.attributes.items():
                        if dataset_name not in root[group]:
                            log("Couldn't find dataset:"
                                " {} in group {}".format(key, dataset_name),
                                cat='error')
                        output[key] = np.array(root[group][dataset_name])
                    output['name'] = path.stem + ('-' + group)
                    outputs.append(output)
            return [self.object(**output) for output in outputs]

        except TypeError as err:
            log(str(err))

        except Exception as err:
            log("Ignoring '{}'"
                " Caught ({})".format(path.name, type(err).__name__),
                cat='warning')
            log(str(err))


def write_hdf5_file(fname, attributes):
    """Writes matrix data to hdf5 file"""
    with h5py.File(fname, 'w') as root:
        for key in attributes.keys():
            data = attributes[key]
            dset = root.create_dataset(key,
                                       data.shape,
                                       data.dtype,
                                       compression='gzip')
            dset[...] = data[...]


def standard_figure(extent=(0.5, 2.5, 0.5, 2.5), figsize=(4, 4), dpi=300):
    """Standard initialization for HS fingerprint figures"""
    plt.style.use('seaborn-white')
    fig = plt.figure(figsize=figsize, dpi=dpi)
    cmap = viridis
    axes = fig.add_subplot(111, xlim=extent[0:2], ylim=extent[2:4])
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
    return fig, cmap, axes, extent


def hexbin_plotfile(xvals, yvals, fname='outhex.png',
                    kind='log', nbins=100):
    """Create a hexagonally binned HS fingerprint"""
    extent = [min(0.5, np.min(xvals)), max(2.5, np.max(xvals)*1.2),
              min(0.5, np.min(yvals)), max(2.5, np.max(yvals)*1.2)]
    fig, cmap, _, extent = standard_figure(extent=extent)

    plt.hexbin(xvals, yvals, mincnt=1, gridsize=nbins,
               cmap=cmap, bins=kind,
               extent=extent)

    fig.savefig(fname, bbox_inches='tight')
    plt.clf()
    plt.close()


# FILE FUNCTIONS
def plotfile(xvals, yvals, fname='out.png', kind='linear', nbins=100):
    """ Construct a histogram, plot it, then write the image as a png
    to a given filename."""

    # Not sure why but we have linear and log as options
    if kind == 'linear':
        hist, xedges, yedges = bin_data(xvals, yvals, bins=nbins)

    # NEED TO ROTATE AXES
    hist = np.rot90(hist)
    hist = np.flipud(hist)

    fig, cmap, _, _ = standard_figure()

    # this is the key step, the rest is formatting
    plt.pcolormesh(xedges, yedges, hist,
                   cmap=cmap,
                   norm=mpl.colors.LogNorm())

    plt.grid(b=True, which='major', axis='both')

    # SAVE TO PNG WITH LITTLE TO NO BORDER
    fig.savefig(fname, bbox_inches='tight')

    # CLOSE UP
    plt.clf()
    plt.close()


def write_sa_file(fname, surfaces):
    """ Write out the surface area contribution info to a given
    file."""
    with open(fname, 'w') as surface_file:
        for surface in surfaces:
            line = '{0}, {1}'.format(surface.name, surface.formula)
            if not surface.contributions:
                line = line + '--Nil--'
            else:
                for interaction, percent in surface.contributions.items():
                    line = line + ', '
                    line = line + '{} = {:.2%}'.format(interaction,
                                                       percent)
            surface_file.write(line + '\n')


def write_mat_file(fname, mat, names, clusters):
    """ Writes clustering info to file """
    matrix_dict = defaultdict()
    log("Writing clustering data to  '{}'".format(str(fname)))
    matrix_dict['matrix'] = np.array(mat)
    matrix_dict['names'] = np.array(names, dtype='<S64')
    matrix_dict['clusters'] = np.array(clusters)
    write_hdf5_file(fname, matrix_dict)


def batch_process(files, reader, procs=4):
    """
    Takes a list of files and a reader class, and
    constructs instances of the given class from data
    in the files using ProcessPoolExecutor
    """
    nfiles = len(files)
    with Timer() as time:
        vals = []

        with concurrent.futures.ProcessPoolExecutor(procs) as executor:
            return_values = [executor.submit(reader.read, file) for file in files]
            for value in concurrent.futures.as_completed(return_values):
                vals.append(value.result())

        vals = [x for sublist in vals for x in sublist if sublist is not None]

        n_success = len(vals)
        if n_success < nfiles:
            err = 'Skipped {0} files due to errors'.format(nfiles - n_success)
            log(err, cat='error')

    log("Read {} files in {}".format(nfiles, time))

    return sorted(vals, key=lambda x: x.name)
