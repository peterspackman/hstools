import os
from collections import OrderedDict

from .config import log
from .datafile import (
        DataFileReader,
        SurfaceData,
        batch_process,
        write_sa_file
)
surface_keys = {'contributions':'surface_contribution', 'formula':'formula'}

def proc_file_sa(fname, order=False, matrix=False):
    """ Process an input file for use in calculating the
    contribution of element -> element interactions on the
    hirshfeld surface """
    ret = None
    cname = get_basename(fname)
    r = readh5file(fname, ["surface_contribution", "formula"])

    sa = r["surface_contribution"]
    formula = str(r["formula"], 'utf-8')

    if matrix:
        sa += np.triu(np.rollaxis(sa, 1), 1)
        sa = np.triu(sa)
        sa /= np.sum(sa)

        ret = Surface(sa, formula, cname)
    else:
        _, contrib_p = get_contrib(sa, order=order)
        ret = Surface(contrib_p, formula, cname)

    return ret


def process_files(files, output=None):
    reader = DataFileReader(surface_keys,
                            SurfaceData)


    surfaces = batch_process(files, reader)

    # If we are writing to file
    if output:
        write_sa_file(fname, surfaces)

    # Otherwise we are printing to stdout
    else:
        for x in surfaces:
            log('Molecular Formula: {0}'.format(str(x.formula, 'utf-8')))
            if x.contributions is None:
                log(' -- Nil--')

            d = OrderedDict(sorted(x.contributions_dict.items(), key=lambda t: t[1]))
            for k, v in iter(d.items()):
                log('{0}: {1:.2%}'.format(k, v))
