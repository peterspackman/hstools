from scipy.spatial import KDTree, ConvexHull
from scipy.special import sph_harm
import trimesh
import numpy as np
import sbf
from .lebedev import lebedev_grid, integrate_values
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

log = logging.getLogger('decompose')

def spherical_to_cartesian(rtp):
    """
    Given an N by 3 array of (r, theta, phi) spherical coordinates
    return an N by 3 array of Cartesian(x, y, z) coordinates.
    """
    xyz = np.zeros(rtp.shape)

    xyz[:, 0] = rtp[:, 0] * np.sin(rtp[:, 1]) * np.cos(rtp[:, 2])
    xyz[:, 1] = rtp[:, 0] * np.sin(rtp[:, 1]) * np.sin(rtp[:, 2])
    xyz[:, 2] = rtp[:, 0] * np.cos(rtp[:, 1])

    return xyz


def interpolate(idxs, norms):
    return np.mean(norms[idxs])

def mean_radius(pts):
    center = np.mean(pts, axis=1)
    pts = pts.transpose() - center
    d2 = pts[:, 0] **2 + pts[:, 1] **2 + pts[:, 2] **2
    norms = np.sqrt(d2)
    mean_radius = np.mean(norms)
    return mean_radius

def mean_radius2(pts):
    return np.mean(np.linalg.norm(pts.transpose() -
                   np.mean(pts, axis=1), axis=1))


def describe_surface(filename, degree=131, export_mesh=False):
    name = Path(filename).stem
    f = sbf.File(filename)
    f.read()
    pts = f['vertices'].data.transpose()

    log.debug('Loaded vertex data')
    center = np.mean(pts, axis=0)
    # shift to be centered about the origin
    pts = pts - center 
    if export_mesh:
        faces = f['faces'].data.transpose() - 1
        mesh = trimesh.Trimesh(vertices=pts, faces=faces)
        trimesh.io.export.export_mesh(mesh, 'original.ply')

    # this is faster for some reason than np.apply_along_axis
    norms = np.sqrt(pts[:,0] **2 + pts[:, 1] **2 + pts[:, 2] **2)
    mean_radius = np.mean(norms)
    pts = pts / mean_radius
    norms = norms / mean_radius 
    pts_normalized = pts / np.reshape(norms, (pts.shape[0], 1))
    log.debug('Normalized points')
    grid = lebedev_grid(degree=degree)
    grid_cartesian= spherical_to_cartesian(np.c_[np.ones(grid.shape[0]), grid[:, 1], grid[:, 0]])
    if export_mesh:
        faces = ConvexHull(pts_normalized).simplices
        mesh = trimesh.Trimesh(vertices=pts_normalized, faces=faces)
        trimesh.io.export.export_mesh(mesh, 'sphere.ply')
        faces = ConvexHull(grid_cartesian).simplices
        mesh = trimesh.Trimesh(vertices=grid_cartesian, faces=faces)
        trimesh.io.export.export_mesh(mesh, 'sphere-lebedev.ply')

    log.debug('Constructing tree')
    tree = KDTree(pts_normalized)
    log.debug('Done')
    log.debug('Interpolating values')
    nn = tree.query(grid_cartesian, 1)
    log.debug('Done')
    values = np.array([interpolate(idxs, norms) for idxs in nn[1]])
    return name, sht(values, grid)

def reconstruct_surface(coeffs, l_max=20, degree=131):
    """Reconstruct the HS by distorting a sphere"""
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
    coefficients = np.zeros((l_max + 1)*(l_max + 1), dtype=np.complex128)
    lm = 0
    for l in range(0, l_max + 1):
        for m in range(-l, l + 1):
            vals = sph_harm(m, l, grid[:, 0], grid[:, 1]) * values
            coefficients[lm] = integrate_values(grid, vals)
            lm += 1
    return coefficients

def export_mesh(verts, faces, filename='output.ply'):
    mesh = trimesh.Trimesh(vertices=verts, faces=faces)
    trimesh.io.export.export_mesh(mesh, filename)

def read_radius(arrfile, root):
    parent = arrfile.parent.name
    sbf_file = Path(root, parent, arrfile.with_suffix('.sbf').name)
    f = sbf.File(sbf_file)
    f.read()
    radius = mean_radius(f['vertices'].data)
    return radius



def main():
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument('directory')
    parser.add_argument('-l', '--lmax', default=20, type=int,
                        help='Maximum angular momentum')
    parser.add_argument('--log-file', default=None,
                        help='Log to file instead of stdout')
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
    
