"""
Command line interface for writing information about SA contribution
based on interatomic interactions on the surface
"""
from .config import log
from .datafile import (
    DataFileReader,
    SurfaceData,
    batch_process,
)

SURFACE_KEYS = {'contributions': 'surface_contribution'} 


def process_files(files, output=None):
    """
    Given a list of hdf5 files (Path objects), read the atomic
    contributions to the Hirshfeld surface area, and report the results.

    Returns a list of SurfaceData objects
    """
    reader = DataFileReader(SURFACE_KEYS,
                            SurfaceData)

    surfaces = batch_process(files, reader)

    # TODO
    if output:
        pass

    for surface in surfaces:
        log('Molecular Formula: {0}'.format(surface.name.split('-')[-1]))
        if surface.contributions is None:
            log(' -- Nil--')
        else:
            for key, value in surface.contributions_dict.items():
                log('{0}: {1:.2%}'.format(key, value))

    return surfaces
