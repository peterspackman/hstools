# Core imports
import sys
import os
import glob
from multiprocessing.pool import ThreadPool
import time
# Library imports
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import progressbar as pb
# Local imports
import hist
import calc
import pack.ccalc as ccalc
from data import widgets
from data import log


# HELPER AND WRAPPER FUNCTIONS
def get_vals(lines, t=np.float64, index=0):
    """ return a list, interpreting the nth word on each line
        as the value to be stored based on the data type"""
    return [t(line.split()[index]) for line in lines]


def surface_helper(args):
    """ Helper function for map_async to proc_file_sa"""
    fname, restrict, order = args
    return proc_file_sa(fname, restrict, order=order)


def hist_helper(args):
    """ Helper function for map_async on proc_file_hist"""
    fname, res, save_figs = args
    return proc_file_hist(fname, resolution=res, save_figs=save_figs)


def harmonics_helper(args):
    fname, metric = args
    return proc_file_harmonics(fname, metric=metric)


def readcxsfile_c(fname):
    """ A wrapper around ccalc.readcxsfile """
    try:
        r = ccalc.readcxsfile(fname)
    except ccalc.error, e:
        log('Problem in {0}: {1}'.format(fname, e))
        return None
    di, de, p, harmonics = r
    formula, vertices, indices, internal, external = p
    # Strip the unnecessary quotes and spaces from the line
    formula = formula.split('\"')[1]
    np.require(internal, requirements=['O'])
    return di, de, (formula, vertices, indices, internal, external), harmonics


# FILE FUNCTIONS
def plotfile(x, y, fname='out.png', type='linear', nbins=10):
    """ Construct a histogram, plot it, then write the image as a png
    to a given filename."""

    # Not sure why but we have linear and log as options
    if(type == 'linear'):
        H, xedges, yedges = hist.bin_data(x, y, bins=nbins)
    else:
        H, xedges, yedges = hist.bin_data_log(x, y, bins=nbins)

    # NEED TO ROTATE AXES
    H = np.rot90(H)
    H = np.flipud(H)

    fig = plt.figure()

    # this is the key step, the rest is formatting
    plt.pcolormesh(xedges, yedges, H, norm=mpl.colors.LogNorm())

    # set x and y limits, the values to draw ticks, and turn on the grid
    plt.xlim([0.5, 2.5])
    plt.ylim([0.5, 2.5])
    plt.xticks(np.arange(0.5, 2.5, 0.2))
    plt.yticks(np.arange(0.5, 2.5, 0.2))
    plt.grid()

    ticklabelpad = mpl.rcParams['xtick.major.pad']
    # Label our graph axes in nice places!
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


def proc_file_sa(fname, restrict, order=False):
    """ Process an input file for use in calculating the
    contribution of element -> element interactions on the
    hirshfeld surface """
    if not os.path.isfile(fname):
        err = 'Could not open {0} for reading, check to see if file exists'
        log(err.format(fname))
        sys.exit(1)
    r = readcxsfile_c(fname)
    if not r:
        return None
    x, y, a, harmonics = r

    formula, vertices, indices, internal, external = a
    contrib, contrib_p = calc.get_contrib_percentage(vertices, indices,
                                                     internal, external,
                                                     x + y, dp=1,
                                                     restrict=restrict,
                                                     order=order)
    cname = os.path.basename(os.path.splitext(fname)[0])
    return cname, formula, contrib_p


def proc_file_harmonics(fname, metric='dnorm'):
    """ Read a file from fname, collecting the coefficients/invariants
        and returning them as arrays
    """
    if not os.path.isfile(fname):
        err = 'Could not open {0} for reading, check to see if file exists'
        if os.path.isdir(fname):
            err = '{0} appears to be a directory, use --batch'
        log(err.format(fname))
        sys.exit(1)
    r = readcxsfile_c(fname)
    cname = os.path.basename(os.path.splitext(fname)[0])
    if not r:
        return None
    _, _, _, harmonics = r
    return (harmonics, cname)


def proc_file_hist(fname, resolution=10, save_figs=False):
    """ Read a file from fname, generate a histogram and potentially write
        the png of it to file.
        FUTURE: offer the restrictions of internal and external atoms
        using internal[] and external[] along with indices[]
    """
    if not os.path.isfile(fname):
        err = 'Could not open {0} for reading, check to see if file exists'
        if os.path.isdir(fname):
            err = '{0} appears to be a directory, use --batch'
        log(err.format(fname))
        sys.exit(1)
    r = readcxsfile_c(fname)
    if not r:
        return None
    x, y, a, _ = r

    cname = os.path.basename(os.path.splitext(fname)[0])
    h = hist.bin_data(x, y, resolution)
    if(save_figs):
        outfile = os.path.splitext(fname)[0] + '{0}bins.png'.format(resolution)
        plotfile(x, y, fname=outfile, nbins=resolution)

    return h, cname


def write_sa_file(fname, cnames, formulae, contribs):
    """ Write out the surface area contribution info to a given
    file."""
    with open(fname, 'w') as f:
        for i in range(len(formulae)):
            cname = cnames[i]
            formula = formulae[i]
            contrib_p = contribs[i]
            line = '{0}, {1}'.format(cname, formula)
            if not contrib_p:
                line = line + '--Nil--'
            else:
                for key in sorted(contrib_p, key=lambda key: contrib_p[key]):
                    line = line + ', '
                    line = line + '{0} = {1:.2%}'.format(key, contrib_p[key])
            f.write(line + '\n')


def write_mat_file(fname, mat):
    np.savetxt(fname, mat, fmt="%.4e", delimiter=' ')


# BATCH FUNCTIONS
def batch_hist(dirname, suffix='.cxs', resolution=10,
               save_figs=False, procs=4):
    """Generate n histograms from a directory, returning a list of them
       and their corresponding substance names """
    if not os.path.isdir(dirname):
        err = '{0} does not appear to be a directory'
        log(err.format(dirname))
        sys.exit(1)
    files = sorted(glob.glob(os.path.join(dirname, '*'+suffix)))
    nfiles = len(files)
    if nfiles < 1:
        log('No files to read in {0}'.format(dirname))
        sys.exit(1)
    args = [(fname, resolution, save_figs) for fname in files]

    histograms = []
    names = []
    vals = []
    pbar = pb.ProgressBar(widgets=widgets, maxval=nfiles)
    start_time = time.time()
    pbar.start()
    p = ThreadPool(procs)
    r = p.map_async(hist_helper, args, callback=vals.extend)
    p.close()
    # Doing something I ought not to do, using private members
    # of MapResult to check how many are done. (bad)
    while True:
        if r.ready():
            break
        pbar.update()
        time.sleep(0.2)
    p.join()
    pbar.finish()
    # Strip none values
    vals = [x for x in vals if x is not None]
    if len(vals) < nfiles:
        log('Skipped {0} files due to errors'.format(nfiles - len(vals)))
        nfiles = len(vals)

    if nfiles > 0:
        # unzip the output
        histograms, names = zip(*vals)
        output = 'Reading {0} files took {1:.2} seconds using {2} threads.'
        log(output.format(nfiles, time.time() - start_time, procs))

        return (histograms, names)
    else:
        log('Errors reading all files, exiting.')
        sys.exit(1)


def batch_harmonics(dirname, metric='d_norm', suffix='.cxs', procs=4):
    if not os.path.isdir(dirname):
        err = '{0} does not appear to be a directory'
        log(err.format(dirname))
        sys.exit(1)
    files = sorted(glob.glob(os.path.join(dirname, '*'+suffix)))
    nfiles = len(files)
    if nfiles < 1:
        log('No files to read in {0}'.format(dirname))
        sys.exit(1)
    args = [(fname, metric) for fname in files]

    values = []
    names = []
    vals = []
    # Boilerplate
    pbar = pb.ProgressBar(widgets=widgets, maxval=nfiles)
    start_time = time.time()
    pbar.start()
    p = ThreadPool(procs)
    r = p.map_async(harmonics_helper, args, callback=vals.extend)
    p.close()
    while True:
        if r.ready():
            break
        pbar.update()
        time.sleep(0.2)
    p.join()
    pbar.finish()
    # Strip none values
    vals = [x for x in vals if x is not None]
    if len(vals) < nfiles:
        log('Skipped {0} files due to errors'.format(nfiles - len(vals)))
        nfiles = len(vals)

    if nfiles > 0:
        # unzip the output
        values, names = zip(*vals)
        output = 'Reading {0} files took {1:.2} seconds using {2} threads.'
        log(output.format(nfiles, time.time() - start_time, procs))

        return (values, names)
    else:
        log('Errors reading all files, exiting.')
        sys.exit(1)


def batch_surface(dirname, restrict, suffix='.cxs', procs=4, order=False):
    """ Traverse a directory calculating the surface area contribution"""
    if not os.path.isdir(dirname):
        err = '{0} does not appear to be a directory'
        log(err.format(dirname))
        sys.exit(1)
    files = sorted(glob.glob(os.path.join(dirname, '*'+suffix)))
    nfiles = len(files)
    if nfiles < 1:
        log('No files to read in {0}'.format(dirname))
        sys.exit(1)
    args = [(fname, restrict, order) for fname in files]

    formulae = []
    contribs = []
    vals = []
    pbar = pb.ProgressBar(widgets=widgets, maxval=nfiles)
    start_time = time.time()
    pbar.start()

    p = ThreadPool(procs)
    r = p.map_async(surface_helper, args, callback=vals.extend)
    p.close()

    while True:
        if r.ready():
            break
        pbar.update()
        time.sleep(0.2)
    p.join()
    pbar.finish()
    # Strip none values
    vals = [x for x in vals if x is not None]
    if len(vals) < nfiles:
        log('Skipped {0} files due to errors'.format(nfiles - len(vals)))
        nfiles = len(vals)

    if nfiles > 0:
        cnames, formulae, contribs = zip(*vals)
        output = 'Reading {0} files took {1:.2} seconds using {2} threads.'
        log(output.format(nfiles, time.time() - start_time, procs))
        return (cnames, formulae, contribs)
    else:
        log('Errors reading all files, exiting.')
        sys.exit(1)
