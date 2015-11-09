from plyfile import PlyElement, PlyData
from scipy.special import sph_harm
from collections import namedtuple
import numpy as np
import matplotlib as mpl
from pathlib import Path

from .datafile import DataFileReader
from .cmap import viridis, rviridis

hirshfeld_defaults = {'vertices':'vertices', 'indices':'indices', 'property':'d_norm'}
coefficient_defaults = {'coefficients':'coefficients', 'property_coefficients':'dnorm_coefficients'}

RGB = namedtuple('RGB', 'red green blue')

HirshfeldData = namedtuple('SurfaceData', "name vertices indices property")
CoefficientsData = namedtuple('CoefficientsData', "name coefficients property_coefficients")


def construct_surface(coeff, vol, lmax, radius=1.0):
    """
    Construct the surface from the provided coefficients
    
    """
    r = bounds(coeff, lmax)
    M, N, P  = vol
    # Separation in voxels
    sepx = r*2/M
    sepy = r*2/N
    sepz = r*2/P
    sep = (sepx, sepy, sepz)
    center =  np.array([sepx * M/2, sepy *N/2, sepz *P/2])
    iso_level = 0.0
    volume = inside(coeff, vol, sep, lmax)
    print('Starting marching cubes')
    verts, faces = measure.marching_cubes(volume, iso_level, sep)
    print('Done\nStarting face correction')
    corrected_faces = measure.correct_mesh_orientation(volume, verts, faces, sep)
    print('Done, center: {}'.format(center))
    verts -= center
    return verts, corrected_faces, center


def get_reconstructed_surface(data, lmax, vol=(100,100,100)):
    coeff = data.coefficients
    color_coeff = data.property_coefficients
    print('Marching cubes on {}, using {} volume.'.format(data.name, vol))
    verts, faces, center = construct_surface(data.coefficients, vol, lmax)
    print('Done')
    rtp = cartesianToSpherical(verts)
    color_vals = np.zeros(len(rtp))
    lm = 0
    for l in range(0,lmax+1):
        for m in range(-l,l+1):
            color_vals[:] += (data.property_coefficients[lm] * 
                             sph_harm(m, l, rtp[:,2], rtp[:,1])).real
            lm += 1
    color_vals *= 4*np.pi
    return verts, faces, vertex_colors(color_vals)


def bounds(coeff, lmax):
    t, p = np.linspace(0, 2*np.pi,50), np.linspace(0,np.pi, 50)
    tp = cartesian([t,p])
    r = np.zeros(len(tp))
    r[:] += coeff[0].real \
            * sph_harm(0, 1, tp[:,0],
                       tp[:,1]).real
    val = np.max(r)
    return val * 1.5 if val > 1.0 else 1.5

def cartesian(arrays, out=None):
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in range(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

def inside(coeff, vol, sep, lmax):
    r = 0
    x, y, z = np.indices(vol)
    x = (x - vol[0]/2)*sep[0]
    y = (y - vol[1]/2)*sep[1]
    z = (z - vol[2]/2)*sep[2]

    voxel_data = np.zeros(vol)
    xy =  x[:,:,:]**2 + y[:,:,:]**2

    r = np.sqrt(xy + z[:,:,:]**2)
    theta = np.arctan2(z[:,:,:], np.sqrt(xy)) + np.pi/2
    phi = np.arctan2(y[:,:,:], x[:,:,:]) + np.pi
    constructed = np.zeros(r.shape)

    lm = 0
    print('Heavy loop')
    for l in range(0,lmax+1):
        for m in range(-l,l+1):
            constructed[:,:,:] += (coeff[lm] \
                     * sph_harm(m, l, phi[:,:,:],
                                theta[:,:,:])).real
            lm += 1
    print('finished')
    constructed = r - constructed
    return constructed

def cartesianToSpherical(xyz):
    ptsnew = np.zeros(xyz.shape)
    xy =  xyz[:,0]**2 + xyz[:,1]**2
    ptsnew[:,0] = np.sqrt(xy + xyz[:,2]**2)
    ptsnew[:,1] = np.arctan2(xyz[:,2], np.sqrt(xy)) + np.pi/2
    ptsnew[:,2] = np.arctan2(xyz[:,1], xyz[:,0]) + np.pi
    return ptsnew

def sphericalToCartesian(rtp):
    ptsnew = np.zeros(rtp.shape)
    ptsnew[:,0] = rtp[:,0] * np.sin(rtp[:,1]) * np.cos(rtp[:,2])
    ptsnew[:,1] = rtp[:,0] * np.sin(rtp[:,1]) * np.sin(rtp[:,2])
    ptsnew[:,2] = rtp[:,0] * np.cos(rtp[:,1])
    return ptsnew


def color_function(val, min, max, property='d_norm'):
    startColor = RGB(255,0,0)
    midColor = RGB(255,255,255)
    endColor = RGB(0,0,255)
    LIMIT = 0.0001
    val = (val)/(max -min)
    if val < 0.0:
        factor = 1.0 - val / min
        color = startColor
    else:
        factor = 1.0 - val/ max
        color = endColor

    if factor > 0.0:
        return RGB(int(color.red + (midColor.red - color.red) * factor),
                   int(color.green + (midColor.green - color.green) * factor),
                   int(color.blue + (midColor.blue - color.blue) * factor))
    else:
        return color


def map_viridis(data, norm=None):
    m = mpl.cm.ScalarMappable(norm=norm, cmap=rviridis)
    return m.to_rgba(data)


def vertex_colors(property_values):
    minima = np.min(property_values)
    maxima = np.max(property_values)
    norm = mpl.colors.Normalize(vmin=minima, vmax=maxima)
    rgb = np.ndarray.astype(255*map_viridis(property_values, norm=norm), 'u1')
    return rgb

def get_HS_data(data, property='d_norm'):
    verts = data.vertices
    colors = vertex_colors(data.property)
    faces = data.indices - 1
    return verts, faces, colors


def write_ply_file(verts, faces, colors, output_file='dump.ply'):
    print("Writing file to {}".format(output_file))
    vertices = np.zeros(len(verts),
                        dtype=([('x', 'f4'), ('y', 'f4'), ('z', 'f4'),
                                ('red', 'u1'), ('green', 'u1'),
                                ('blue', 'u1')]))

    vertices[:]['x'] = verts[:,0]
    vertices[:]['y'] = verts[:,1]
    vertices[:]['z'] = verts[:,2]

    vertices[:]['red'] = colors[:,0]
    vertices[:]['green'] = colors[:,1]
    vertices[:]['blue'] = colors[:,2]

    indices = np.zeros(len(faces),
                       dtype=[('vertex_indices', 'i4', (3,))])

    indices['vertex_indices'] = faces[:,:]

    surface_data = PlyData(
        [
            PlyElement.describe(vertices, 'vertex',
                                comments=['surface vertices']),
            PlyElement.describe(indices, 'face')

        ]
    )
    surface_data.write(output_file)


def process_files(files, reconstruct=False, output=None):
    lmax = 10
    for f in files:
        if reconstruct:
            reader = DataFileReader(coefficient_defaults, CoefficientsData)
            data = reader.read(f)
            verts, faces, colors = get_reconstructed_surface(data, lmax)
            output_file = f.stem + '-reconstructed.ply'
        else:
            reader = DataFileReader(hirshfeld_defaults, HirshfeldData)
            data = reader.read(f)
            verts, faces, colors = get_HS_data(data, property='d_norm')
            output_file = f.stem + '-hirshfeld.ply'
        write_ply_file(verts, faces, colors, output_file=output_file)
