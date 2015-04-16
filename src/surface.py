"""
Usage:
    sarlacc surface [options] [<filepattern>...]

Options:

    -i=ATOM, --internal-atom=ATOM  Restrict the closest internal atom
                                   in the histogram
    -e=ATOM, --external-atom=ATOM  Restrict the closest external atom
                                   in the histogram
    --order-important              When classifying surface area,
                                   indicate that H -> O is different
                                   to O -> H. (i.e. order is important)
    -o=FILE, --output=FILE         Write the result to a given file.
                                   Works for both surface and hist modes.
                                   Will save the distance matrix in hist mode,
                                   and write S.A. infor in surface mode

"""
# Core imports
import os
from collections import OrderedDict

# Library imports
from docopt import docopt

# Local imports
from . import calc
from .data import log, log_error, logClosestPair, logFarthestPair
from .fileio import proc_file_sa, batch_surface, write_sa_file


def process_file_list(files, args, procs):
    if len(files) < 1:
        log_error("No files to process.")
        return
    surfaces = batch_surface(files,
               procs=procs,
               order=args['--order-important'])

    # If we are writing to file
    if args['--output']:
        fname = args['--output']
        write_sa_file(fname, surfaces)

    # Otherwise we are printing to stdout
    else:
        for x in surfaces:
            log('Molecular Formula: {0}'.format(x.formula))
            if not x.contributions:
                log(' -- Nil--')

            d = OrderedDict(sorted(x.contributions.items(), key=lambda t: t[1]))
            for k, v in iter(d.items()):
                log('{0}: {1:.2%}'.format(k, v))


def surface_main(argv, procs=4):
    args = docopt(__doc__, argv=argv)

    order = args['--order-important']
    if len(args['<filepattern>']) < 2:
        file_pattern = args['<filepattern>'][0]

        if os.path.isfile(file_pattern):
            # Generate the percentage contribution of each element
            cname, formula, contrib_p = proc_file_sa(file_pattern,
                                                     order=order)
            log('{0} {1}'.format(cname, formula))

            d = OrderedDict(sorted(contrib_p.items(), key=lambda t: t[1]))
            for k, v in iter(d.items()):
                log('{0}: {1:.2%}'.format(k, v))

        elif os.path.isdir(file_pattern):
            from .fileio import glob_directory
            with glob_directory(file_pattern, '*.hdf5') as files:
                process_file_list(files, args, procs)
        else:
            from glob import glob
            files = glob(file_pattern)
            if len(files) < 1:
                log_error("pattern {} matched no files.".format(file_pattern))
            else:
                process_file_list(files, args, procs)

    else:
        process_file_list(args['<filepattern>'], args, procs)
