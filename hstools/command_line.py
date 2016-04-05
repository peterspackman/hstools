"""
Entry point for hstools script
"""
import argparse
import sys
from pathlib import Path
from hstools.config import Timer, log

class CLIParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: {}\n'.format(message))
        self.print_help()
        sys.exit(2)

def read_paths(patterns, suffix, process_files, **kwargs):
    """
    Given a list of patterns and a suffix, resolve all
    paths from the patterns with the matching suffix, then
    call process_files on the resulting list of files with
    kwargs as the arguments.

    """
    files = set()
    needle = '*.{}'.format(suffix)
    for pattern in patterns:
        path = Path(pattern).resolve()
        if path.is_dir():
            matches = list(path.glob(needle))
            if len(matches) < 1:
                log("No files in '{}' matching '{}'".format(
                    path.absolute(),
                    needle),
                    cat='warning')
            files = files.union(matches)
        elif path.is_file():
            files.add(path)

    if len(files) > 0:
        process_files(files, **kwargs)
    else:
        log('No files to process.', cat='warning')


def harmonics(args):
    """
    Cluster based on rotationally invariant spherical harmonic shape
    descriptors.

    This looks for all HDF5 files in the paths given
    and will analyse the invariants, using them as a vector for comparison
    in a hierarchical clustering.
    """
    from hstools.harmonics import process_files
    read_paths(args.paths, args.suffix, process_files,
               no_radius=args.no_radius,
               output=args.output)


def fingerprint(args):
    """
    Cluster based on Hirshfeld fingerprints as descriptors.

    This will look for all HDF5 files in the paths given in order
    to construct hstools fingerprints. The histogram based on the
    d_e/d_i values in these data files will then be used as the vector
    for a hierarchical clustering.
    """
    from hstools.fingerprint import process_files
    read_paths(args.paths, args.suffix, process_files,
               png=args.make_figures, output=args.output)


def surface(args):
    """
    Calculate Hirshfeld surface composition descriptors.

    """
    from hstools.surface import process_files
    read_paths(args.paths, args.suffix, process_files)


def mesh(args):
    """
    Export/generate .ply files from HDF5 data
    """
    from hstools.mesh import process_files
    read_paths(args.paths, args.suffix, process_files,
               reconstruct=args.reconstruct, cmap=args.colormap,
               surface_property=args.property, lmax=args.lmax, output=args.output)


def describe(args):
    """
    TODO
    process CIF files, describing the hirshfeld surface
    """
    from hstools.describe import process_files
    read_paths(args.paths, args.suffix, process_files,
               surface_property=args.property, lmax=args.lmax)


def main():
    """
    Entry point for the program
    """
    parser = CLIParser()
    subparsers = parser.add_subparsers()

    # HARMONICS
    harmonics_p = subparsers.add_parser('harmonics')
    harmonics_p.add_argument('paths', nargs='*')
    harmonics_p.add_argument('--suffix', default='h5',
                             help='suffix for HDF5 files to look for')
    harmonics_p.add_argument('--output', default=None)
    harmonics_p.add_argument('--no-radius', action='store_true')
    harmonics_p.set_defaults(func=harmonics)

    # FINGERPRINT
    fingerprint_p = subparsers.add_parser('fingerprint')
    fingerprint_p.add_argument('paths', nargs='*')
    fingerprint_p.add_argument('--suffix', default='h5',
                               help='suffix for HDF5 files to look for')
    fingerprint_p.add_argument('--make-figures', action='store_true')
    fingerprint_p.add_argument('--output', default=None)
    fingerprint_p.set_defaults(func=fingerprint)

    # SURFACE
    surface_p = subparsers.add_parser('surface')
    surface_p.add_argument('paths', nargs='*')
    surface_p.add_argument('--suffix', default='h5',
                           help='suffix for HDF5 files to look for')
    surface_p.set_defaults(func=surface)

    # MESH
    mesh_p = subparsers.add_parser('mesh')
    mesh_p.add_argument('paths', nargs='*')
    mesh_p.add_argument('--suffix', default='h5',
                        help='suffix for HDF5 files to look for')
    mesh_p.add_argument('--property', default='d_norm',
                        choices={'d_norm', 'd_i', 'd_e', 'curvature'},
                        help='property to decorate HS')
    mesh_p.add_argument('--reconstruct', action='store_true',
                        help='generate HS from coefficients')
    mesh_p.add_argument('--colormap', default='viridis_r',
                        help='map to color the surface based on property')
    mesh_p.add_argument('--lmax', default=9,
                        help='maximum l-value to use for reconstruction')
    mesh_p.add_argument('--output', default=None)
    mesh_p.set_defaults(func=mesh)

    with Timer() as time:
        args = parser.parse_args()
        args.func(args)
    log('Complete {}'.format(time))
