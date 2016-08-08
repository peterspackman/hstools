"""
Command line functionality for export/reconstructing HS
from HDF5 files
"""
from collections import namedtuple
import trimesh
from tqdm import tqdm
from scipy.special import sph_harm
from scipy.spatial import ConvexHull
import numpy as np
import matplotlib as mpl

from .datafile import DataFileReader
from .lebedev import lebedev_grid

HIRSHFELD_DEFAULTS = {'vertices': 'vertices',
                      'indices': 'indices',
                      'property': 'd_norm'}
COEFFICIENT_DEFAULTS = {'coefficients': 'coefficients',
                        'property_coefficients': 'dnorm_coefficients',
                        'radius': 'radius'}

RGB = namedtuple('RGB', 'red green blue')

HirshfeldData = namedtuple('SurfaceData',
                           'name vertices indices property')
CoefficientsData = namedtuple('CoefficientsData',
                              'name coefficients property_coefficients radius')

log = logging.getLogger(__name__)


def cartesian_to_spherical(xyz):
    """
    Given an N by 3 array of (r, theta, phi) spherical coordinates
    return an N by 3 array of Cartesian(x, y, z) coordinates.
    """
    rtp = np.zeros(xyz.shape)
    x_y = xyz[:, 0]**2 + xyz[:, 1]**2
    rtp[:, 0] = np.sqrt(x_y + xyz[:, 2]**2)
    rtp[:, 1] = np.arctan2(xyz[:, 2], np.sqrt(x_y)) + np.pi/2
    rtp[:, 2] = np.arctan2(xyz[:, 1], xyz[:, 0]) + np.pi
    return rtp


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


def map_colors(data, norm=None, cmap='viridis_r'):
    """Map the viridis colormap to an array of data"""
    return mpl.cm.ScalarMappable(norm=norm, cmap=cmap).to_rgba(data)


def vertex_colors(property_values, cmap='viridis_r'):
    """Given an array of property values, return an array of colors"""
    minima = np.min(property_values)
    maxima = np.max(property_values)
    norm = mpl.colors.Normalize(vmin=minima, vmax=maxima)
    rgb = np.ndarray.astype(255 * map_colors(property_values,
                                              norm=norm,
                                              cmap=cmap), 'u1')
    return rgb


def get_hs_data(data, cmap='viridis_r'):
    """Get the HS data out of a DataReader object"""
    verts = data.vertices
    colors = vertex_colors(data.property, cmap=cmap)
    faces = data.indices - 1
    return verts, faces, colors

def get_delaunay_surface(data, lmax):
    """Reconstruct the HS by distorting a sphere"""
    grid = lebedev_grid(degree=131)
    sphere_rtp = np.zeros(grid.shape)
    sphere_rtp[:, 0] = 1.0
    sphere_rtp[:, 1] = grid[:, 1]
    sphere_rtp[:, 2] = grid[:, 0]
    verts = sphere_rtp.copy()
    sphere = spherical_to_cartesian(sphere_rtp)
    radius = np.zeros(len(verts), dtype=np.complex)
    colors = np.zeros(len(verts), dtype=np.complex)
    l_m = 0
    for l in range(0, lmax + 1):
        for m in range(-l, l + 1):
            ylm = sph_harm(m, l, grid[:, 0], grid[:, 1])
            radius[:] += data.coefficients[l_m] * ylm
            colors[:] += data.property_coefficients[l_m] * ylm
            l_m += 1
    verts[:, 0] = radius[:].real
    verts = spherical_to_cartesian(verts) * data.radius
    colors *= 4 * np.pi
    faces = ConvexHull(sphere).simplices
    return verts, faces, vertex_colors(colors.real)


def process_files(files, reconstruct=False, lmax=9, **kwargs):
    """
    Given a list of HDF5 files, export/reconstruct their
    Hirshfeld surface (with colouring) to a .stl file.
    """
    cmap = kwargs.get('cmap', 'viridis_r')
    if reconstruct:
        log.info('Reconstructing HS from coefficients')
        reader = DataFileReader(COEFFICIENT_DEFAULTS, CoefficientsData)
        suffix = '-reconstructed-HS.stl'
        get_surface = lambda g: get_delaunay_surface(g, lmax)
    else:
        log.info('Extracting HS from vertex/face data')
        HIRSHFELD_DEFAULTS['property'] = kwargs.get('property', 'd_norm')
        reader = DataFileReader(HIRSHFELD_DEFAULTS, HirshfeldData)
        suffix = '-HS.stl'
        get_surface = lambda g: get_hs_data(g, cmap=cmap)

    for path in tqdm(files, unit='file', leave=True):
        data = reader.read(path)
        for group in tqdm(data, unit='HS', leave=True):
            verts, faces, colors = get_surface(group)
            output_file = group.name + suffix
            mesh = trimesh.Trimesh(vertices=verts, faces=faces)
            mesh.visual.vertex_colors = colors
            trimesh.io.export.export_mesh(mesh, output_file)
