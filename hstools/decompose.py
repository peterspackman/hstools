from scipy.spatial import cKDTree as KDTree, ConvexHull
from scipy.special import sph_harm
import numpy as np
import sbf
from .lebedev import lebedev_grid, integrate_values
import logging
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

log = logging.getLogger('decompose')


def spherical_to_cartesian(rtp):
    """
    Given an N by 3 array of (r, theta, phi) spherical coordinates
    return an N by 3 array of Cartesian(x, y, z) coordinates.

    Arguments:
    rtp -- set of r, theta, phi coordinates
    """
    xyz = np.zeros(rtp.shape)

    xyz[:, 0] = rtp[:, 0] * np.sin(rtp[:, 1]) * np.cos(rtp[:, 2])
    xyz[:, 1] = rtp[:, 0] * np.sin(rtp[:, 1]) * np.sin(rtp[:, 2])
    xyz[:, 2] = rtp[:, 0] * np.cos(rtp[:, 1])

    return xyz


def _interpolate(idxs, norms):
    return np.mean(norms[idxs])


def shift_to_origin(pts):
    """reoriginate a set of points to be centered about [0, 0, 0]

    Arguments:
    pts -- set of points to reoriginate
    """
    center = np.mean(pts, axis=0)
    return pts - center


def mean_radius(pts, reoriginate=False):
    """Calculate the mean radius (distance to origin) of a set
    of vertices

    Arguments:
    pts -- set of points to calculate the mean norm

    Keyword arguments:
    reoriginate -- shift the points to be centered about [0,0,0] first
    (default False)"""
    if reoriginate:
        pts = shift_to_origin(pts)
    d2 = pts[:, 0] ** 2 + pts[:, 1] ** 2 + pts[:, 2] ** 2
    norms = np.sqrt(d2)
    mean_radius = np.mean(norms)
    return mean_radius


def values_from_grid(vals, ix):
    res = np.array([_interpolate(idxs, vals) for idxs in ix])
    return res


def describe_surface(filename, degree=131, properties=None):
    """Given an SBF, describe the set of vertices inside
    using spherical harmonics. Will scale the mesh to be of unit
    mean radius.

    Arguments:
    filename -- name of the SBF to open

    Keyword arguments:
    degree -- the degree of lebedev grid to use for integration
    (default 131)
    (default False)
    """
    if not properties:
        properties = []

    name = Path(filename).stem
    f = sbf.File(filename)
    f.read()
    pts = f['vertices'].data.transpose()

    log.debug('Loaded vertex data')
    center = np.mean(pts, axis=0)
    # shift to be centered about the origin
    pts -= center

    # this is faster for some reason than np.apply_along_axis
    norms = np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2 + pts[:, 2] ** 2)
    mean_radius = np.mean(norms)
    pts /= mean_radius
    norms /= mean_radius
    pts_normalized = pts / np.reshape(norms, (pts.shape[0], 1))
    log.debug('Normalized points')
    grid = lebedev_grid(degree=degree)
    grid_cartesian = spherical_to_cartesian(
            np.c_[np.ones(grid.shape[0]), grid[:, 1], grid[:, 0]])
    log.debug('Constructing tree')
    tree = KDTree(pts_normalized)
    log.debug('Done')
    log.debug('Interpolating values')
    nn = tree.query(grid_cartesian, 1)
    log.debug('Done')
    values = values_from_grid(norms, nn[1])
    property_arrays = {'shape': values}

    for prop in properties:
        property_arrays[prop] =  values_from_grid(f[prop].data, nn[1])
    property_coefficients = { name: sht(value, grid) for name, value in property_arrays.items()}

    return name, mean_radius, property_coefficients


def combination_desc(filename, degree=131):
    """Given an SBF, describe the set of vertices and their esp using sht.
    Will scale the mesh to be of unit mean radius.

    Arguments:
    filename -- name of the SBF to open

    Keyword arguments:
    degree -- the degree of lebedev grid to use for integration
    (default 131)
    (default False)
    """
    name = Path(filename).stem
    f = sbf.File(filename)
    f.read()
    pts = f['vertices'].data.transpose()

    log.debug('Loaded vertex data')
    center = np.mean(pts, axis=0)
    # shift to be centered about the origin
    pts -= center

    # this is faster for some reason than np.apply_along_axis
    norms = np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2 + pts[:, 2] ** 2)
    mean_radius = np.mean(norms)
    pts /= mean_radius
    norms /= mean_radius
    pts_normalized = pts / np.reshape(norms, (pts.shape[0], 1))
    log.debug('Normalized points')
    grid = lebedev_grid(degree=degree)
    grid_cartesian = spherical_to_cartesian(
            np.c_[np.ones(grid.shape[0]), grid[:, 1], grid[:, 0]])
    log.debug('Constructing tree')
    tree = KDTree(pts_normalized)
    log.debug('Done')
    log.debug('Interpolating values')
    nn = tree.query(grid_cartesian, 1)
    log.debug('Done')
    shape = values_from_grid(norms, nn[1])
    esp =  values_from_grid(f['electric_potential'].data, nn[1])
    esp_max = np.max(np.abs(esp))
    esp /= esp_max
    combined = np.zeros(esp.shape, dtype=np.complex128)
    combined = shape + np.conj(esp)
    coefficients = sht(combined, grid)
    return name, mean_radius, esp_max, coefficients



def reconstruct_surface(coeffs, l_max=20, degree=131):
    """Reconstruct the HS by distorting a spherical mesh, generated
    from a lebedev grid.

    Arguments:
    coeffs -- the set of spherical harmonic coefficients

    Keyword arguments:
    l_max -- maximum angular momenta to reconstruct to (default 20)
    degree -- grid degree to use (see lebedev grids)
    """
    grid = lebedev_grid(degree=degree)
    # grid[:, 0] goes from 0 -> 2 PI
    # grid[:, 1] goes from 0 -> PI
    rtp = np.c_[np.ones(grid.shape[0]), grid[:, 1], grid[:, 0]]
    verts = rtp.copy()
    sphere = spherical_to_cartesian(rtp)
    radius = np.zeros(len(verts), dtype=np.complex)
    l_m = 0
    for l in range(0, l_max + 1):
        for m in range(-l, l + 1):
            ylm = sph_harm(m, l, grid[:, 0], grid[:, 1])
            radius[:] += coeffs[l_m] * ylm
            l_m += 1
    verts[:, 0] = radius[:].real
    verts = spherical_to_cartesian(verts)
    faces = ConvexHull(sphere).simplices
    return verts, faces


def sht(values, grid, l_max=20):
    """Perform a spherical harmonic transform given a grid and a set of values

    Arguments:
    values -- set of scalar function values associated with grid points
    grid -- set of grid points

    Keyword arguments:
    l_max -- maximum angular momenta to integrate to (default 20)
    """
    coefficients = np.zeros((l_max + 1)*(l_max + 1), dtype=np.complex128)
    lm = 0
    for l in range(0, l_max + 1):
        for m in range(-l, l + 1):
            vals = sph_harm(m, l, grid[:, 0], grid[:, 1]) * values
            coefficients[lm] = integrate_values(grid, vals)
            lm += 1
    return coefficients


def main():
    """Read through all sbf files in a directory, writing
    numpy arrays of sht coefficients for each shape encountered.
    """
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument('directory')
    parser.add_argument('-l', '--lmax', default=20, type=int,
                        help='Maximum angular momentum')
    parser.add_argument('--log-file', default=None,
                        help='Log to file instead of stdout')
    parser.add_argument('--suffix', '-s', default='.sbf',
                        help='File suffix to find sbf files')
    parser.add_argument('--log-level', default='INFO',
                        help='Log level')
    parser.add_argument('--jobs', '-j', default=4, type=int,
                        help='Number of parallel jobs to run')
    parser.add_argument('--output-directory', '-o', default='.',
                        help='Directory to store output numpy arrays')
    args = parser.parse_args()
    if args.log_file:
        logging.basicConfig(filename=args.log_file, level=args.log_level)
    else:
        logging.basicConfig(level=args.log_level)
    log.info('Starting %s, output: %s', args.directory, args.output_directory)
    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)

    paths = list(Path(args.directory).glob('*.sbf'))
    log.info('%d paths to process', len(paths))
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(describe_surface, path) for path in paths]
        for f in as_completed(futures):
            name, coeffs = f.result()
            array_file = os.path.join(args.output_directory, name)
            log.info('Saving array to %s', array_file)
            np.save(array_file, coeffs)
    log.info('Finished %s', args.directory)

if __name__ == '__main__':
    main()
