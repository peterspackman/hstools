# Core imports
import sys
import os
import glob
import time
import string
# Library imports
from joblib import Parallel, delayed
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
# Local imports
import hist
import visual
import pack.cio as cio
# A temporary variable for the formatting of the histogram plots
ticklabelpad = mpl.rcParams['xtick.major.pad']


def get_vals(lines, t=np.float64, index=0):
    """ return a numpy array, interpreting the first word on each line
        as the value to be stored """
    return [t(line.split()[index]) for line in lines]


def readcxsfile_c(fname):
    di, de, p = cio.readcxsfile(fname)
    formula, vertices, indices, internal, external = p
    # Strip the unnecessary quotes and spaces from the line
    formula = formula.split('\"')[1]
    return di, de, (formula, vertices, indices, internal, external)


def plotfile(x, y, fname='out.png', type='linear', nbins=10):
    """ Construct a histogram, plot it, then write the image as a png
    to a given filename.
    """

    # Not sure why i've bothered with linear and log as options
    if(type == 'linear'):
        H, xedges, yedges = hist.bin_data(x, y, bins=nbins)
    else:
        H, xedges, yedges = hist.bin_data_log(x, y, bins=nbins)

    # NEED TO ROTATE AXES
    H = np.rot90(H)
    H = np.flipud(H)

    Hmasked = np.ma.masked_where(H == 0, H)

    fig = plt.figure()

    # this is the key step, the rest is formatting
    plt.pcolormesh(xedges, yedges, H, norm=mpl.colors.LogNorm())

    # set x and y limits, the values to draw ticks, and turn on the grid
    plt.xlim([0.5, 2.5])
    plt.ylim([0.5, 2.5])
    plt.xticks(np.arange(0.5, 2.5, 0.2))
    plt.yticks(np.arange(0.5, 2.5, 0.2))
    plt.grid()

    # Label our graph axes in nice places !
    plt.annotate(r'$d_i$', fontsize=20, xy=(1, 0), xytext=(5, - ticklabelpad),
                 ha='left', va='top', xycoords='axes fraction',
                 textcoords='offset points')
    plt.annotate(r'$d_e$', fontsize=20, xy=(0, 1.02),
                 xytext=(5, - ticklabelpad), ha='right',
                 va='bottom', xycoords='axes fraction',
                 textcoords='offset points')

    # SAVE TO PNG WITH LITTLE TO NO BORDER
    fig.savefig(fname, bbox_inches='tight')

    # CLOSE UP
    plt.clf()
    plt.close()


def process_file(fname, resolution=10, write_png=False, i=None, e=None):
    """ Read a file from fname, generate a histogram and potentially write
        the png of it to file. i restricts internal atom, e restricts external
    """
    if not os.path.isfile(fname):
        err = 'Could not open {0} for reading, check to see if file exists'
        print err.format(fname)
        sys.exit(1)

    x, y, a = readcxsfile_c(fname)

    cname = os.path.basename(os.path.splitext(fname)[0])
    h = hist.bin_data(x, y, resolution)
    if(write_png):
        outfile = os.path.splitext(fname)[0] + '{0}bins.png'.format(resolution)
        plotfile(x, y, fname=outfile, nbins=resolution)

    return h, cname


def batch_process(dirname, suffix='.cxs', resolution=10,
                  write_png=False, threads=4, i=None, e=None):
    """Generate n histograms from a directory, returning a list of them
       and their corresponding substance names
       Note that the 'threads' here are actually processes"""
    if not os.path.isdir(dirname):
        err = '{0} does not appear to be a directory'
        print err.format(dirname)
        sys.exit(1)
    files = glob.glob(os.path.join(dirname, '*'+suffix))

    histograms = []
    names = []

    start_time = time.time()

    # Now this is some ugly indentation and syntax!
    vals = Parallel(n_jobs=threads,
                    verbose=3)(delayed(process_file)(f,
                                                     resolution=resolution,
                                                     write_png=write_png,
                                                     i=i, e=e)
                               for f in files)

    # unzip the output
    if (i or e):
        vals = [(h, n) for h, n in vals if h and n]
    histograms, names = zip(*vals)
    output = 'Reading {0} files took {1:.2} seconds using {2} threads.'
    print output.format(len(files),time.time() - start_time, threads)

    return (histograms, names)
