# Core imports
import concurrent.futures
import glob
import os
import sys
import time
# Library imports
from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
import progressbar as pb
# Local imports
from data import log
from data import widgets
import calc
import hist

nmdims = {"vertices": 3, "indices": 3, "coefficients": 2}
ndtypes = {"indices": np.int32, "atoms_inside_surface": np.int32,
           "atoms_outside_surface": np.int32,
           "d_e_face_atoms": np.int32, "d_i_face_atoms": np.int32,
           "unit_cell": np.dtype((str, 3))}
numerical = {"unit_cell": True}


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


def split_text(s):
    from itertools import groupby
    for k, g in groupby(s, str.isalpha):
        yield ''.join(list(g))


# Read a cxs file getting values specified in names
# returns a dictionary contianing the values
def readcxsfile(fname, attributes):
    outputs = {}
    reading = None
    count = 0
    dtype = np.float64
    expectedVals = 1
    with open(fname) as f:
        for line in f:

            if reading and count > 0:
                arr = outputs[reading]
                try:
                    x = np.fromstring(line, dtype=dtype,
                                      count=expectedVals, sep=' ')
                    if(x.size > 1):
                        if(len(arr.shape)) < 2:
                            log('Retrieved more values than expected????')
                            raise ValueError
                        arr[arr.shape[0] - count, :] = x
                    else:
                        arr[arr.size - count] = x
                except:
                    x = next(split_text(line))
                    arr[arr.size - count] = x
                count -= 1

            elif line.startswith('begin '):
                name = line.split()[1]

                if name in attributes:
                    reading = name
                    count = int(line.split()[2])

                    if name in ndtypes:
                        dtype = ndtypes[name]
                    else:
                        dtype = np.float64

                    if name in nmdims:
                        expectedVals = nmdims[name]
                        outputs[name] = np.zeros((count,
                                                  expectedVals),
                                                 dtype=dtype)

                    else:
                        expectedVals = 1
                        outputs[name] = np.zeros(count, dtype=dtype)
            elif line.startswith("   formula = "):
                outputs["formula"] = line.split('"')[1]

    return outputs


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
    r = readcxsfile(fname, ["vertices", "indices", "atoms_inside_surface",
                            "atoms_outside_surface", "d_i_face_atoms",
                            "d_e_face_atoms", "unit_cell", "d_e", "d_i"])
    formula = r["formula"]
    cname = os.path.basename(os.path.splitext(fname)[0])

    x = r["vertices"]
    y = r["indices"]
    ai = r["atoms_inside_surface"]
    ao = r["atoms_outside_surface"]
    di = r["d_i_face_atoms"]
    de = r["d_e_face_atoms"]
    uc = r["unit_cell"]
    distances = r["d_e"] + r["d_i"]
    external = uc[ao[de - 1] - 1]
    internal = uc[ai[di - 1] - 1]

    _, contrib_p = calc.get_contrib_percentage(x, y, internal,
                                               external, distances,
                                               restrict=restrict, order=order)
    return (cname, formula, contrib_p)


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
    r = readcxsfile(fname, ["coefficients", "invariants"])
    cname = os.path.basename(os.path.splitext(fname)[0])
    if not r:
        return None
    harmonics = (r["coefficients"], r["invariants"])
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
    r = readcxsfile(fname, ["d_e", "d_i"])
    if not r:
        return None

    x = r["d_i"]
    y = r["d_e"]

    cname = os.path.basename(os.path.splitext(fname)[0])
    h = hist.bin_data(x, y, resolution)
    if(save_figs):
        outfile = os.path.splitext(fname)[0] + '{0}bins.png'.format(resolution)
        plotfile(x, y, fname=outfile, nbins=resolution)

    ret = (h, cname)
    return ret


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

    with concurrent.futures.ProcessPoolExecutor(procs) as executor:
        fs = [executor.submit(hist_helper, arg) for arg in args]
        for i, f in enumerate(concurrent.futures.as_completed(fs)):
            pbar.update(i)
            a = f.result()
            vals.append(a)

    pbar.finish()
    # Strip none values
    vals = [x for x in vals if x is not None]
    if len(vals) < nfiles:
        log('Skipped {0} files due to errors'.format(nfiles - len(vals)))
        nfiles = len(vals)

    vals = sorted(vals, key=lambda val: val[1])

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

    with concurrent.futures.ProcessPoolExecutor(procs) as executor:
        fs = [executor.submit(harmonics_helper, arg) for arg in args]
        for i, f in enumerate(concurrent.futures.as_completed(fs)):
            pbar.update(i)
            vals.append(f.result())

    pbar.finish()
    # Strip none values
    vals = [x for x in vals if x is not None]
    if len(vals) < nfiles:
        log('Skipped {0} files due to errors'.format(nfiles - len(vals)))
        nfiles = len(vals)

    vals = sorted(vals)

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

    with concurrent.futures.ProcessPoolExecutor(procs) as executor:
        fs = [executor.submit(surface_helper, arg) for arg in args]
        for i, f in enumerate(concurrent.futures.as_completed(fs)):
            pbar.update(i)
            vals.append(f.result())

    pbar.finish()

    # Strip none values
    vals = [x for x in vals if x is not None]
    if len(vals) < nfiles:
        log('Skipped {0} files due to errors'.format(nfiles - len(vals)))
        nfiles = len(vals)

    vals = sorted(vals)

    if nfiles > 0:
        cnames, formulae, contribs = zip(*vals)
        output = 'Reading {0} files took {1:.2} seconds using {2} threads.'
        log(output.format(nfiles, time.time() - start_time, procs))
        return (cnames, formulae, contribs)
    else:
        log('Errors reading all files, exiting.')
        sys.exit(1)
