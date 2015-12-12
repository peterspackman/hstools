import click
import os.path
import matplotlib
from pathlib import Path

from libsarlacc.harmonics import modes as harmonics_modes
from libsarlacc.config import Timer, log

# Library imports

# Local imports
available_commands = ['harmonics', 'fingerprint',
                      'surface', 'mesh']


def read_paths(patterns, suffix, process_files, **kwargs):
    """
    Given a list of patterns and a suffix, resolve all
    paths from the patterns with the matching suffix, then
    call process_files on the resulting list of files with
    kwargs as the arguments.

    """
    files = []

    match = '*.{}'.format(suffix)
    for pattern in patterns:
        p = Path(pattern)
        if p.is_dir():
            matches = list(p.glob(match))
            if len(matches) <= 0:
                log("No HDF5 files in '{}'".format(p.absolute()),
                    cat='warning')
            files += matches
        elif p.is_file():
            files.append(p)

    if len(files) > 0:
        process_files(files, **kwargs)
    else:
        log('No files to process.', cat='error')


@click.group()
@click.option('--version/-v', default=False, is_flag=True)
def cli(version):
    if version:
        log('VERSION INFO')
        sys.exit(0)


@cli.command()
@click.option('--output', default=None, is_flag=False,
              help='write out clustering to given HDF5 file name')
@click.option('--property', type=click.Choice(list(harmonics_modes.keys())),
              default='shape',
              help='property on which to base the clustering')
@click.option('--suffix', default='h5',
              help='suffix for HDF5 files to look for')
@click.option('--no-radius', default=False,
              help="don't use mean radius in the clustering")
@click.argument('paths', type=click.Path(exists=True), nargs=-1)
def harmonics(paths, output, property, suffix, no_radius):
    """
    Cluster based on rotationally invariant spherical harmonic shape
    descriptors.

    This looks for all HDF5 files in the paths given
    and will analyse the invariants, using them as a vector for comparison
    in a hierarchical clustering.
    """
    from libsarlacc.harmonics import process_files
    with Timer() as t:
        read_paths(paths, suffix, process_files,
                   no_radius=no_radius, mode=property, output=output)
    log('Complete {}'.format(t))


@cli.command()
@click.option('--png / -p', default=False,
              help='output all Hirshfeld fingerprints to PNG image files.')
@click.option('--suffix', default='h5',
              help='suffix for HDF5 files to look for')
@click.option('--output', default='fingerprints.h5', is_flag=False,
              help='write out clustering to given HDF5 file name')
@click.argument('paths', nargs=-1)
def fingerprint(paths, png, suffix, output):
    """
    Cluster based on Hirshfeld fingerprints as descriptors.

    This will look for all HDF5 files in the paths given in order
    to construct Hirshfeld fingerprints. The histogram based on the
    d_e/d_i values in these data files will then be used as the vector
    for a hierarchical clustering.
    """
    from libsarlacc.fingerprint import process_files
    with Timer() as t:
        read_paths(paths, suffix, process_files,
                   png=png, output=output)
        log('Complete {}'.format(t))


@cli.command()
@click.option('--suffix', default='h5',
              help='suffix for HDF5 files to look for')
@click.argument('paths', nargs=-1)
def surface(paths, suffix):
    """
    Calculate Hirshfeld surface composition descriptors.

    """
    from libsarlacc.surface import process_files
    with Timer() as t:
        read_paths(paths, suffix, process_files)
    log('Complete {}'.format(t))


@cli.command()
@click.option('--reconstruct / -',
              help='generate the surface from the coefficients')
@click.option('--suffix', default='h5',
              help='suffix for HDF5 files to look for')
@click.option('--colors', default='viridis_r',
              type=click.Choice(matplotlib.pyplot.colormaps()),
              help='colormap for property on the surface')
@click.option('--property', default='d_norm',
              type=click.Choice(['d_norm', 'd_i', 'd_e', 'curvature']))
@click.argument('paths', nargs=-1)
def mesh(paths, suffix, reconstruct, colors, property):
    """
    Export/generate .ply files from HDF5 data
    """
    from libsarlacc.mesh import process_files
    with Timer() as t:
        read_paths(paths, suffix, process_files,
                   reconstruct=reconstruct, cmap=colors,
                   property=property)
    log('Complete {}'.format(t))
