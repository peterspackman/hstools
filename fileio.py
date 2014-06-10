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
# A temporary variable for the formatting of the histogram plots
ticklabelpad = mpl.rcParams['xtick.major.pad']


def get_vals(lines, t=np.float64):
    """ return a numpy array, interpreting the first word on each line
        as the value to be stored """
    return [t(line.split()[0]) for line in lines]


def readcxsfile(fname):
    """ Hacky way to find the de_vals and di_vals  in a cxs file
        Ideally I'd write a file parser especially for cxs files
        but this is not yet done
    """
    devals = []
    vertices = []
    indices = []
    external = []
    internal = []
    formula = ""

    with open(fname) as f:
        get_count = lambda x: int(x.split()[2])
        get_vertices = lambda x: [float(y) for y in x.split()]
        get_indices = lambda x: [int(y) for y in x.split()]
        content = f.readlines()

        for n in range(len(content)):
            # Basically going through line by line matching
            # patterns and storing data contained within
            line = content[n]

            # FORMULA
            if line.startswith('   formula'):
                formula = content[n].split('\"')[1]

            # VERTICES
            if line.startswith('begin vertices'):
                r = get_count(line)
                vertices = np.zeros((r, 3))
                for i in range(r):
                    n = n + 1
                    vertices[i] = get_vertices(content[n])

            # INDICES
            if line.startswith('begin indices'):
                r = get_count(line)
                indices = np.zeros((r, 3), dtype=np.int32)
                for i in range(r):
                    n = n + 1
                    indices[i] = get_indices(content[n])

            # D_E VALUES
            if line.startswith('begin d_e '):
                r = get_count(line)
                n = n + 1
                devals = np.array(get_vals(content[n:n+r]))
                n = n + r

            # D_I VALUES
            if line.startswith('begin d_i '):
                r = get_count(line)
                n = n + 1
                divals = np.array(get_vals(content[n:n+r]))
                n = n + r

            # D_I FACE ATOM SYMBOLS
            if line.startswith('begin d_i_face_chemical'):
                r = get_count(line)
                n = n + 1
                internal = np.array(get_vals(content[n:n+r], t=str))
                n = n + r

            # D_E ATOM SYMBOLS
            if line.startswith('begin d_e_face_chemical'):
                r = get_count(line)
                n = n + 1
                external = np.array(get_vals(content[n:n+r], t=str))
                n = n + r

    fail1 = type(devals) is list or type(divals) is list
    fail2 = type(internal) is list or type(external) is list
    fail3 = type(vertices) is list or type(indices) is list
    # If we have a problem. i.e. de_vals or di_vals will be empty
    if fail1 or fail2 or fail3:
        print 'FATAL: missing values'
        print 'Input file is likely missing necessary data from tonto'
        sys.exit(0)
    return divals, devals, (formula, vertices, indices, internal, external)


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
    # in essence, the x-axis is di_values (internal) while y is d_e
    x, y, a = readcxsfile(fname)
    formula, _, _, internal, external = a
    if i:
        for ind in range(internal.size):
            if not (internal[ind] == i):
                x[ind] = 0.
                y[ind] = 0.
    if e:
        for ind in range(external.size):
            if not (external[ind] == e):
                x[ind] = 0.
                y[ind] = 0.
    if i or e:
        x = x[np.nonzero(x)]
        y = y[np.nonzero(y)]

    cname = os.path.basename(os.path.splitext(fname)[0])
    if x.size < 1 or y.size < 1:
        print '{0} has no {1} -> {2} interactions'.format(fname, i, e)
        return
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
    output = 'Reading {0} files took {1} seconds using {2} threads.'
    print output.format(len(files), time.time() - start_time, threads)

    return (histograms, names)
