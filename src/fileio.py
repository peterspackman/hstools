# Core imports
from collections import defaultdict, namedtuple
from contextlib import contextmanager
import concurrent.futures
import glob
import os
import sys
import time

# Library imports
from matplotlib import pyplot as plt
import h5py
import matplotlib as mpl
import numpy as np
import progressbar as pb
import seaborn as sns

# Local imports
from . import calc
from . import data
from .data import log, log_traceback, logger

Harmonic = namedtuple('Harmonic', 'coefficients invariants name')
Surface = namedtuple('Surface', 'contributions formula name')
Hist = namedtuple('Hist', 'histogram name')


class EmptyDirException(Exception):
    pass


class FilesSkippedException(Exception):
    pass


def dict_vals(d, *keys):
    values = ()
    for key in keys:
        if key not in d:
            values += (None,)
        else:
            values += (d[key],)
    return values


@contextmanager
def glob_directory(d, pattern):
    try:
        yield sorted(glob.glob(os.path.join(d, pattern)))
    except EmptyDirException as e:
        logger.error(e)
        sys.exit(1)
    except Exception as e:
        logger.exception(e)


def readh5file(fname, attributes):
    outputs = defaultdict()
    with h5py.File(fname, 'r') as f:
        for a in attributes:
            if a not in f:
                logger.error("Couldn't find dset: {} in {}".format(a, f))
            outputs[a] = f[a].value
    return outputs


def load_saved_data(fname):
    attributes = ['matrix', 'names']
    data = readh5file(fname, attributes)
    return data['matrix'], data['names']


# Takes a dictionary of dataset_name
def writeh5file(fname, attributes):
    with h5py.File(fname, 'w') as f:
        for k in attributes.keys():
            data = attributes[k]
            dset = f.create_dataset(k,
                                    data.shape,
                                    data.dtype,
                                    compression='gzip')
            dset[...] = data[...]


def dir_file_join(d, f):
    return [d, f].join('/')


def standard_figure(figsize=(9, 9), dpi=400):
    f = plt.figure(figsize=figsize, dpi=dpi)
    sns.set(style='white')
    cmap = mpl.cm.jet
    extent = [0.5, 2.5, 0.5, 2.5]
    ax = f.add_subplot(111, xlim=extent[0:2], ylim=extent[2:4])
    plt.xticks(np.arange(0.5, 2.5, 0.2))
    plt.yticks(np.arange(0.5, 2.5, 0.2))
    plt.grid(b=True, which='major', axis='both')

    ticklabelpad = mpl.rcParams['xtick.major.pad']

    # Label our graph axes in nice places!
    plt.annotate(r'$d_i$', fontsize=20, xy=(1, 0), xytext=(5, - ticklabelpad),
                 ha='left', va='top', xycoords='axes fraction',
                 textcoords='offset points')
    plt.annotate(r'$d_e$', fontsize=20, xy=(0, 1.02),
                 xytext=(5, - ticklabelpad), ha='right',
                 va='bottom', xycoords='axes fraction',
                 textcoords='offset points')
    return f, cmap, ax, extent


def kde_plotfile(x, y, fname='outhex.png', type='linear'):

    f, cmap, ax, extent = standard_figure()

    sns.kdeplot(x, y, cmap=cmap, shade=True, ax=ax)
    ax.collections[0].set_alpha(0)

    f.savefig(fname, bbox_inches='tight')
    plt.clf()
    plt.close()


def hexbin_plotfile(x, y, fname='outhex.png', type='linear', nbins=10):
    f, cmap, ax, extent = standard_figure()

    image = plt.hexbin(x, y, mincnt=1, gridsize=nbins,
                       cmap=cmap, bins='log',
                       extent=extent)

    f.savefig(fname, bbox_inches='tight')
    plt.clf()
    plt.close()


# FILE FUNCTIONS
def plotfile(x, y, fname='out.png', type='linear', nbins=10):
    """ Construct a histogram, plot it, then write the image as a png
    to a given filename."""

    # Not sure why but we have linear and log as options
    if(type == 'linear'):
        H, xedges, yedges = calc.bin_data(x, y, bins=nbins)
    else:
        H, xedges, yedges = calc.bin_data_log(x, y, bins=nbins)

    # NEED TO ROTATE AXES
    H = np.rot90(H)
    H = np.flipud(H)

    f, cmap, ax, extent = standard_figure()

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


def get_basename(fname):
    return os.path.basename(os.path.splitext(fname)[0])


def proc_file_sa(fname, order=False):
    """ Process an input file for use in calculating the
    contribution of element -> element interactions on the
    hirshfeld surface """
    ret = None
    try:
        cname = get_basename(fname)
        r = readh5file(fname, ["surface_contribution", "formula"])

        sa = r["surface_contribution"]
        formula = str(r["formula"], 'utf-8')

        contrib, contrib_p = calc.get_contrib(sa, order=order)
        ret = Surface(contrib_p, formula, cname)

    except Exception as e:
        logger.warning('Skipping {0} => {1}'.format(fname, str(e)))
    finally:
        return ret


def proc_file_harmonics(fname, property='shape'):
    """ Read a file from fname, collecting the coefficients/invariants
        and returning them as arrays
    """
    try:
        name = get_basename(fname)
        r = readh5file(fname, ["coefficients", "invariants",
                       "dnorm_coefficients", "dnorm_invariants",
                       "curvature_coefficients", "curvature_invariants"])

        if property == 'curvature':
            return Harmonic(r['curvature_coefficients'],
                            r['curvature_invariants'],
                            name)
        elif property == 'dnorm':
            return Harmonic(r["dnorm_coefficients"],
                            r["dnorm_invariants"],
                            name)
        else:
            return Harmonic(r["coefficients"],
                            r["invariants"],
                            name)
    except Exception as e:
        logger.error('Skipping {0} => {1}'.format(fname, str(e)))


def proc_file_hist(fname, resolution=10, save_figs=False):
    """ Read a file from fname, generate a histogram and potentially write
        the png of it to file.
    """
    ret = None
    try:
        cname = get_basename(fname)
        r = readh5file(fname, ["d_e", "d_i"])
        x, y = r["d_i"], r["d_e"]
        h = calc.bin_data(x, y, resolution)

        if(save_figs):
            prefix = os.path.splitext(fname)[0]
            outfile = prefix + 'bins.png'
            plotfile(x, y, fname=outfile, nbins=resolution)
            hexbin_plotfile(x,
                            y,
                            fname='{}-hex.png'.format(prefix),
                            nbins=resolution)
            kde_plotfile(x, y, fname='{}-kde.png'.format(prefix))

        ret = Hist(h, cname)
    except Exception as e:
        print(e)
        logger.warning('Skipping {}'.format(fname))
    finally:
        return ret


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
                for key in sorted(contributions, key=lambda key: contributions[key]):
                    line = line + ', '
                    line = line + '{0} = {1:.2%}'.format(key, contributions[key])
            f.write(line + '\n')


def write_mat_file(fname, mat, names, clusters):
    d = defaultdict()
    d['matrix'] = mat
    d['names'] = names
    d['clusters'] = clusters
    writeh5file(fname, d)


# BATCH FUNCTIONS

# Takes a list of arguments (args), and a function, and calls
# that function with that list of args using a Processpoolexecutor
def batch_process(args, function, procs=4, msg='Processing: ', progress=True):
    vals = []

    if progress:
        pbar = pb.ProgressBar(widgets=data.getWidgets(msg),
                              maxval=len(args))
        pbar.start()

    with concurrent.futures.ProcessPoolExecutor(procs) as executor:
        fs = [executor.submit(function, *arg) for arg in args]
        for i, f in enumerate(concurrent.futures.as_completed(fs)):
            if progress:
                pbar.update(i)
            vals.append(f.result())

    if progress:
        pbar.finish()

    vals = [x for x in vals if x is not None]
    if len(vals) < len(args):
        err = 'Skipped {0} files due to errors'.format(len(args)- len(vals))
        logger.warning(err)

    return sorted(vals, key=lambda x: x.name)


def batch_hist(files, resolution=10,
               save_figs=False, procs=4):
    """Generate n histograms from a directory, returning a list of them
       and their corresponding substance names """

    args = [(fname, resolution, save_figs) for fname in files]

    return batch_process(args,
                         proc_file_hist,
                         procs=procs,
                         msg='Reading files: ')


def batch_harmonics(files, property='shape' , suffix='.hdf5', procs=4):

    args = [(fname, property) for fname in files]

    return batch_process(args,
                         proc_file_harmonics,
                         procs=procs,
                         msg='Reading files: ')


def batch_surface(files, suffix='.hdf5', procs=4, order=False):
    """ Traverse a directory calculating the surface area contribution"""

    args = [(fname, order) for fname in files]

    return  batch_process(args,
                          proc_file_sa,
                          procs=procs,
                          msg='Processing files: ')
