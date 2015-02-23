"""
Usage:
    sarlacc surface [options] [<filepattern>]

Options:

    -i=ATOM, --internal-atom=ATOM  Restrict the closest internal atom
                                   in the histogram
    -e=ATOM, --external-atom=ATOM  Restrict the closest external atom
                                   in the histogram
    -j=FILE, --json=FILE           Dump the dendrogram tree to JSON

    --no-restrict                  Don't restrict the surface area
                                   values to only those closer than
                                   Van Der Waal's Radii
    --order-important              When classifying surface area,
                                   indicate that H -> O is different
                                   to O -> H. (i.e. order is important)
    -o=FILE, --output=FILE         Write the result to a given file.
                                   Works for both surface and hist modes.
                                   Will save the distance matrix in hist mode,
                                   and write S.A. infor in surface mode

"""
import os
from collections import OrderedDict

from docopt import docopt
from . import calc
from .data import log
from .fileio import proc_file_sa, batch_surface, write_sa_file
from .modes import logClosestPair, logFarthestPair


def surface_main(argv, procs=4):
    args = docopt(__doc__, argv=argv)

    restrict = not args['--no-restrict']
    order = args['--order-important']
    if os.path.isfile(args['<filepattern>']):
        fname = args['<filepattern>']
        # Generate the percentage contribution of each element
        cname, formula, contrib_p = proc_file_sa(fname, restrict,
                                                 order=order)
        log('{0} {1}'.format(cname, formula))

        d = OrderedDict(sorted(contrib_p.items(), key=lambda t: t[1]))
        for k, v in iter(d.items()):
            log('{0}: {1:.2%}'.format(k, v))

    elif os.path.isdir(args['<filepattern>']):
        dirname = args['<filepattern>']
        cnames, formulae, contribs = batch_surface(dirname, restrict,
                                                   procs=procs,
                                                   order=order)
        if restrict:
            log("Restricted interactions using CCDC Van Der Waal's Radii")
        # If we are writing to file
        if args['--output']:
            fname = args['--output']
            write_sa_file(fname, cnames, formulae, contribs)
        # Otherwise we are printing to stdout
        else:
            for i in range(len(formulae)):
                formula = formulae[i]
                contrib_p = contribs[i]
                log('Molecular Formula: {0}'.format(formula))
                if not contrib_p:
                    log(' -- Nil--')

                d = OrderedDict(sorted(contrib_p.items(), key=lambda t: t[1]))
                for k, v in iter(d.items()):
                    log('{0}: {1:.2%}'.format(k, v))
