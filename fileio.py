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


def readcxsfile(fname):
    """ Hacky way to find the de_vals and di_vals  in a cxs file
        Ideally I'd write a file parser especially for cxs files
        but this is not yet done
    """
    i = 0
    devals = []
    divals = []
    atoms = []
    de_face_atoms = []
    di_face_atoms = []
    atoms_outside = []
    atoms_inside = []
    formula = ""

    with open(fname) as f:

        content = f.readlines()

        for n in range(len(content)):
            # Basically going through line by line matching
            # patterns and saving data
            # FORMULA
            if content[n].startswith('   formula'):
                formula = content[n].split('\"')[1]
            # LIST OF ATOMS
            if content[n].startswith('begin unit_cell'):
                words = content[n].split()
                for i in range(int(words[2])):
                    n = n + 1
                    words = content[n].split()
                    atoms.append(words[1])
            # D_E VALUES
            if content[n].startswith('begin d_e '):
                words = content[n].split()
                r = int(words[2])
                # Update the line counter
                x = np.zeros(r)
                for i in range(r):
                    n = n + 1
                    x[i] = float(content[n])
                devals = x
            # D_I VALUES
            if content[n].startswith('begin d_i '):
                words = content[n].split()
                r = int(words[2])
                # Update the line counter
                x = np.zeros(r)
                for i in range(r):
                    n = n + 1
                    x[i] = float(content[n])
                divals = x
            # D_I FACE ATOMS
            if content[n].startswith('begin d_i_face_atoms'):
                words = content[n].split()
                r = int(words[2])
                x = np.zeros(r, dtype=np.int32)
                for i in range(r):
                    n = n + 1
                    # NOT SURE HOW THEY'RE INDEXED IN THE FILE
                    x[i] = int(content[n])
                di_face_atoms = x
            # D_E FACE ATOMS
            if content[n].startswith('begin d_e_face_atoms'):
                words = content[n].split()
                r = int(words[2])
                x = np.zeros(r, dtype=np.int32)
                for i in range(r):
                    n = n + 1
                    x[i] = int(content[n])
                de_face_atoms = x
            if content[n].startswith('begin atoms_outside'):
                words = content[n].split()
                r = int(words[2])
                x = np.zeros(r, dtype=np.int32)
                for i in range(r):
                    n = n + 1
                    words = content[n].split()
                    x[i] = int(words[0])
                atoms_outside = x
            if content[n].startswith('begin atoms_inside'):
                words = content[n].split()
                r = int(words[2])
                x = np.zeros(r, dtype=np.int32)
                for i in range(r):
                    n = n + 1
                    words = content[n].split()
                    x[i] = int(words[0])
                atoms_inside = x

    l = de_face_atoms.size
    external = np.chararray(l, itemsize=2)  # Array of element names (2chars)
    internal = np.chararray(l, itemsize=2)
    for i in range(l):
        external[i] = atoms[atoms_outside[de_face_atoms[i] - 1] - 1]
        internal[i] = atoms[atoms_inside[di_face_atoms[i] - 1] - 1]

    # We have a problem. i.e. de_vals or di_vals will be empty
    if not devals.any() or not divals.any():
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


def process_file(fname, resolution=10, write_png=False):
    """ Read a file from fname, generate a histogram and potentially write
        the png of it to file
    """
    # note that a is unused currently due to recent change to readcxsfile
    # in essence, the x-axis is di_values (internal) while y is d_e
    x, y, a = readcxsfile(fname)
    cname = os.path.basename(os.path.splitext(fname)[0])
    h = hist.bin_data(x, y, resolution)

    if(write_png):
        outfile = os.path.splitext(fname)[0] + '{0}bins.png'.format(resolution)
        plotfile(x, y, fname=outfile, nbins=resolution)

    return h, cname


def batch_process(dirname, suffix='.cxs', resolution=10,
                  write_png=False, threads=4):
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
                                                     write_png=write_png)
                               for f in files)

    # unzip the output
    histograms, names = zip(*vals)

    output = 'Reading {0} files took {1} seconds using {2} threads.'
    print output.format(len(files), time.time() - start_time, threads)

    return (histograms, names)
