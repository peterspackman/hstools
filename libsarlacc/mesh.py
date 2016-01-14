from plyfile import PlyElement, PlyData
from scipy.special import sph_harm
from collections import namedtuple
import numpy as np
import matplotlib as mpl
import skimage.measure as measure
from scipy.spatial import ConvexHull
from pathlib import Path
from numba import jit

from .datafile import DataFileReader
from .config import log
from .lebedev import lebedev_grid

hirshfeld_defaults = {'vertices': 'vertices',
                      'indices': 'indices',
                      'property': 'd_norm'}
coefficient_defaults = {'coefficients': 'coefficients',
                        'property_coefficients': 'dnorm_coefficients',
                        'radius': 'radius'}

RGB = namedtuple('RGB', 'red green blue')

HirshfeldData = namedtuple('SurfaceData',
                           'name vertices indices property')
CoefficientsData = namedtuple('CoefficientsData',
                              'name coefficients property_coefficients radius')



def cartesian_to_spherical(xyz):
    """
    Given an N by 3 array of (r, theta, phi) spherical coordinates
    return an N by 3 array of Cartesian(x, y, z) coordinates.
    """
    rtp = np.zeros(xyz.shape)
    xy = xyz[:, 0]**2 + xyz[:, 1]**2
    rtp[:, 0] = np.sqrt(xy + xyz[:, 2]**2)
    rtp[:, 1] = np.arctan2(xyz[:, 2], np.sqrt(xy)) + np.pi/2
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


def color_function(val, min_value, max_value, surface_property='d_norm'):
    startColor = RGB(255, 0, 0)
    midColor = RGB(255, 255, 255)
    endColor = RGB(0, 0, 255)
    LIMIT = 0.0001
    val = (val)/(max_value - min_value)
    if val < 0.0:
        factor = 1.0 - val / min_value
        color = startColor
    else:
        factor = 1.0 - val / max_value
        color = endColor

    if factor > 0.0:
        return RGB(int(color.red + (midColor.red - color.red) * factor),
                   int(color.green + (midColor.green - color.green) * factor),
                   int(color.blue + (midColor.blue - color.blue) * factor))
    else:
        return color


def map_viridis(data, norm=None, cmap='viridis_r'):
    m = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    return m.to_rgba(data)


def vertex_colors(property_values, cmap='viridis_r'):
    minima = np.min(property_values)
    maxima = np.max(property_values)
    norm = mpl.colors.Normalize(vmin=minima, vmax=maxima)
    rgb = np.ndarray.astype(255 * map_viridis(property_values,
                                              norm=norm,
                                              cmap=cmap), 'u1')
    return rgb


def get_HS_data(data, cmap='viridis_r'):
    verts = data.vertices
    colors = vertex_colors(data.property, cmap=cmap)
    faces = data.indices - 1
    return verts, faces, colors


def write_ply_file(verts, faces, colors, output_file='dump.ply'):
    log("Writing file to {}".format(output_file))
    vertices = np.zeros(len(verts),
                        dtype=([('x', 'f4'), ('y', 'f4'), ('z', 'f4'),
                                ('red', 'u1'), ('green', 'u1'),
                                ('blue', 'u1')]))

    vertices[:]['x'] = verts[:, 0]
    vertices[:]['y'] = verts[:, 1]
    vertices[:]['z'] = verts[:, 2]

    vertices[:]['red'] = colors[:, 0]
    vertices[:]['green'] = colors[:, 1]
    vertices[:]['blue'] = colors[:, 2]

    indices = np.zeros(len(faces),
                       dtype=[('vertex_indices', 'i4', (3,))])

    indices['vertex_indices'] = faces[:, :]

    surface_data = PlyData(
        [
            PlyElement.describe(vertices, 'vertex',
                                comments=['surface vertices']),
            PlyElement.describe(indices, 'face')

        ]
    )
    surface_data.write(output_file)


def get_delaunay_surface(data, lmax):
    grid = lebedev_grid(degree=131)
    sphere_rtp = np.zeros(grid.shape)
    sphere_rtp[:, 0] = 1.0
    sphere_rtp[:, 1] = grid[:, 1]
    sphere_rtp[:, 2] = grid[:, 0]
    verts = sphere_rtp.copy()
    sphere = spherical_to_cartesian(sphere_rtp)
    r = np.zeros(len(verts), dtype=np.complex)
    colors = np.zeros(len(verts), dtype=np.complex)
    lm = 0
    for l in range(0, lmax+1):
        for m in range(-l, l+1):
            ylm = sph_harm(m, l, grid[:, 0], grid[:, 1])
            r[:] += data.coefficients[lm] * ylm
            colors[:] += data.property_coefficients[lm] * ylm
            lm += 1
    verts[:, 0] = r[:].real
    verts = spherical_to_cartesian(verts) * data.radius
    colors *= 4*np.pi
    faces = ConvexHull(sphere).simplices
    return verts, faces, vertex_colors(colors.real)


def process_files(files, reconstruct=False, output=None,
                  property='d_norm', cmap='viridis_r', lmax=9):
    """
    Given a list of HDF5 files, export/reconstruct their
    Hirshfeld surface (with colouring) to a .ply file.
    """
    output_file = output
    for f in files:
        if reconstruct:
            reader = DataFileReader(coefficient_defaults, CoefficientsData)
            data = reader.read(f)
            verts, faces, colors = get_delaunay_surface(data, lmax)
            if not output:
                output_file = f.stem + '-reconstructed.ply'
        else:
            hirshfeld_defaults['property'] = property
            reader = DataFileReader(hirshfeld_defaults, HirshfeldData)
            data = reader.read(f)
            verts, faces, colors = get_HS_data(data,
                                               cmap=cmap)
            if not output:
                output_file = f.stem + '-hirshfeld.ply'
        write_ply_file(verts, faces, colors, output_file=output_file)
