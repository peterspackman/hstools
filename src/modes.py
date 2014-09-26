# Core imports
import sys
import time
# Library imports
import numpy as np
# Local imports
from .data import log
from . import calc
from . import fileio as fio

test_f = {'sp': calc.spearman_roc,
          'kt': calc.kendall_tau,
          'hd': calc.hdistance}
test_names = {'sp': 'Spearman rank order coefficient',
              'kt': "Kendall's Tau",
              'hd': 'naive histogram distance',
              'dv': 'Euclidean distance between invariants'}

def logClosestPair(mat, names):
    np.fill_diagonal(mat, np.inf)
    x = np.nanargmin(mat)
    minind = (x//len(names), x % len(names))
    log('Closest pair: {0}, d= {1}'.format((names[minind[0]],
                                            names[minind[1]]),
                                            mat[minind]))
    np.fill_diagonal(mat, 0.0)


def hist_main(args):
    mtest = test_f[args['--test']]
    tname = test_names[args['--test']]
    start_time = time.time()
    procs = int(args['--procs'])

    bins = int(args['--bins'])
    save_figs = args['--save-figures']

    if args['<file>']:
        fname = args['<file>']
        if not save_figs:
            log('Not saving figure, so this \
                        command will have no output')
        h, name = fio.proc_file_hist(fname, resolution=bins,
                                     save_figs=save_figs)

    elif args['<dir>']:
        dirname = args['<dir>']
        dendrogram = args['--dendrogram']
        method = args['--method']
        distance = float(args['--distance'])
        # Program is being run to batch process a directory of cxs files
        histograms, names = fio.batch_hist(dirname, resolution=bins,
                                           save_figs=save_figs,
                                           procs=procs)

        log('Generating matrix using {0}'.format(tname))
        mat = calc.get_dist_mat(histograms, test=mtest, threads=procs*2)
        if args['--output']:
            fname = args['--output']
            fio.write_mat_file(fname, mat)

        logClosestPair(mat, names)

        calc.cluster(mat, names, tname, dump=args['--json'],
                     dendrogram=dendrogram,
                     method=method,
                     distance=distance)
    footer(start_time)


def harmonics_main(args):
    mtest = calc.dvalue
    tname = test_names['dv']
    start_time = time.time()
    procs = int(args['--procs'])

    if args['<file>']:
        fname = args['<file>']
        if not fname.endswith('.cxs'):
            log('WARNING: {0} does not have .cxs extension'.format(fname))
        values, cname = fio.proc_file_harmonics(fname)
        coefficients, invariants = values
        log(cname)
        log(invariants)

    if args['<dir>']:
        dendrogram = args['--dendrogram']
        method = args['--method']
        distance = float(args['--distance'])
        dirname = args['<dir>']
        values, names = fio.batch_harmonics(dirname, procs=procs)
        log('Generating matrix using: "{0}"'.format(tname))
        coefficients, invariants = zip(*values)
        mat = calc.get_dist_mat(invariants, test=mtest, threads=procs*2)
        if args['--output']:
            fname = args['--output']
            fio.write_mat_file(fname, mat)
        logClosestPair(mat, names)



        calc.cluster(mat, names, tname, dendrogram=dendrogram,
                     method=method, distance=distance)


    footer(start_time)


def surface_main(args):
    start_time = time.time()
    procs = int(args['--procs'])

    restrict = not args['--no-restrict']
    order = args['--order-important']
    if args['<file>']:
        fname = args['<file>']
        if not fname.endswith('.cxs'):
            log('WARNING: {0} does not have .cxs extension'.format(fname))
        # Generate the percentage contribution of each element
        cname, formula, contrib_p = fio.proc_file_sa(fname, restrict,
                                                     order=order)
        log('{0} {1}'.format(cname, formula))

        for key in sorted(contrib_p, key=lambda key: contrib_p[key]):
            log('{0}: {1:.2%}'.format(key, contrib_p[key]))

    elif args['<dir>']:
        dirname = args['<dir>']
        cnames, formulae, contribs = fio.batch_surface(dirname, restrict,
                                                       procs=procs,
                                                       order=order)
        if restrict:
            log("Restricted interactions using CCDC Van Der Waal's Radii")
        # If we are writing to file
        if args['--output']:
            fname = args['--output']
            fio.write_sa_file(fname, cnames, formulae, contribs)
        # Otherwise we are printing to stdout
        else:
            for i in range(len(formulae)):
                formula = formulae[i]
                contrib_p = contribs[i]
                log('Molecular Formula: {0}'.format(formula))
                if not contrib_p:
                    log(' -- Nil--')
                for key in sorted(contrib_p,
                                  key=lambda key: contrib_p[key]):
                    log('{0}: {1:.2%}'.format(key,
                                              contrib_p[key]))

    footer(start_time)


def footer(start_time):
    log('Process complete! Took {0:.2} s'.format(time.time() - start_time))
    sys.exit(0)
