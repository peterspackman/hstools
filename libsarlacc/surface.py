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

            d = OrderedDict(sorted(x.contributions_dict.items(), key=lambda t: t[0]))
            for k, v in iter(d.items()):
                log('{0}: {1:.2%}'.format(k, v))
