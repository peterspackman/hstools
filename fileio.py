# Core imports
import sys
import os
import glob
import time
from StringIO import StringIO
# Library imports
from joblib import Parallel, delayed
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import from_levels_and_colors
# Local imports
import hist

# A temporary variable for the formatting of the histogram plots
ticklabelpad = mpl.rcParams['xtick.major.pad']


def get_vals(lines, t=float):
    """ return a numpy array, interpreting the first word on each line
        as the value to be stored """
    return  [t(line.strip()[0]) for line in lines]


def readcxsfile(fname):
    """ Hacky way to find the de_vals and di_vals  in a cxs file
        Ideally I'd write a file parser especially for cxs files
        but this is not yet done
    """
    devals = []
    divals = []
    atoms = []
    de_face_atoms = []
    di_face_atoms = []
    atoms_outside = []
    atoms_inside = []
    formula = ""

    with open(fname) as f:
        get_count = lambda x: int(x.split()[2])
        content = f.readlines()

        for n in range(len(content)):
            # Basically going through line by line matching
            # patterns and storing data contained within
            line = content[n]

            # FORMULA
            if line.startswith('   formula'):
                formula = content[n].split('\"')[1]

            # LIST OF ATOMS
            if line.startswith('begin unit_cell'):
                r = get_count(line)
                n = n + 1
                atoms = get_vals(content[n:n+r],t=str)
                n = n + r

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

            # D_I FACE ATOMS
            if line.startswith('begin d_i_face_atoms'):
                r = get_count(line)
                n = n + 1
                di_face_atoms = np.array(get_vals(content[n:n+r],t=int))
                n = n + r

            # D_E FACE ATOMS
            if line.startswith('begin d_e_face_atoms'):
                r = get_count(line)
                n = n + 1
                de_face_atoms = np.array(get_vals(content[n:n+r],t=int))
                n = n + r

            # ATOMS OUTSIDE SURFACE
            if line.startswith('begin atoms_outside'):
                r = get_count(line)
                n = n + 1
                atoms_outside = np.array(get_vals(content[n:n+r],t=int))
                n = n + r

            # ATOMS INSIDE SURFACE
            if line.startswith('begin atoms_inside'):
                r = get_count(line)
                n = n + 1
                atoms_inside = np.array(get_vals(content[n:n+r],t=int))
                n = n + r


    l = de_face_atoms.size  # NUMBER OF POINTS TO DEAL WITH
    external = np.chararray(l, itemsize=2)  # Array of element names (2chars)
    internal = np.chararray(l, itemsize=2)

    # Here we resolve the (nested) references to atomic symbols
    for i in range(l):
        # The indices have to be decremented due to fortran indexing
        external[i] = atoms[atoms_outside[de_face_atoms[i] - 1] - 1]
        internal[i] = atoms[atoms_inside[di_face_atoms[i] - 1] - 1]
        # NOTE something awry with the indexing (wrong di/de mapped)

    # We have a problem. i.e. de_vals or di_vals will be empty
    if not (devals.any() and divals.any()):
        print 'FATAL: missing either d_e or d_i values'
        print 'Input file is likely missing necessary data from tonto'
        sys.exit(0)

    return divals, devals, (formula, internal, external)


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
    oldx = x.size
    oldy = y.size
    formula, internal, external = a
    if i:
        for ind in range(x.size):
            if not (internal[ind] == i):
                x[ind] = 0.
                y[ind] = 0.
    if e:
        for ind in range(y.size):
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
    vals = [i for i in vals if i]
    histograms, names = zip(*vals)
    output = 'Reading {0} files took {1} seconds using {2} threads.'
    print output.format(len(files), time.time() - start_time, threads)

    return (histograms, names)
