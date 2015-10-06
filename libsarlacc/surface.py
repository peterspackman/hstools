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
from .data import log, log_error
from .fileio import batch_surface, write_sa_file


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
    files = []

    for file_pattern in args['<filepattern>']:
        if os.path.isfile(file_pattern):
            files.append(file_pattern)

        elif os.path.isdir(file_pattern):
            from .fileio import glob_directory
            new_files = glob_directory(file_pattern, '*.h5')
            files += new_files
            if len(new_files) < 1:
                log_error("{} no files matching {}.".format(file_pattern, '*.h5'))
        else:
            from glob import glob
            new_files = glob(file_pattern)
            files += new_files
            if len(new_files) < 1:
                log_error("pattern {} matched no files.".format(file_pattern))

    if len(files) > 0:
        process_file_list(files, args, procs)
    else:
        log_error('No files found...')
