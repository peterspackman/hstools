#! /usr/bin/env python3
import os
import sys
from cffi import FFI
import numpy as np
import argparse
from tqdm import tqdm

ffi = FFI()
ffi.cdef('''
// variable setting methods
void tonto_set_basis_directory(const char *input_string, const int length);
void tonto_read_cif(const char *input_string, const int length);
void tonto_close_cif();

//always call these
void tonto_initialize();
void tonto_cleanup();
double res;
int l_max, n_molecules;

//do methods
void tonto_hs_info(const int i, int *n_vertices, int *n_faces);
void tonto_hs_arrays(double *, int *, int, int);
void tonto_print_state();
''')

if sys.platform == 'darwin':
    suffix = 'dylib'
else:
    suffx = 'so'

tonto_location = '/Users/prs/repos/hs/build/src'
tonto = ffi.dlopen('{}/libtontoapi.{}'.format(tonto_location, suffix))

def initialize(args):
    tonto.tonto_set_basis_directory(args.basis.encode("utf-8"),
                                len(args.basis))
    tonto.l_max = args.lmax
    tonto.res = args.resolution
    tonto.tonto_initialize()

def read_cif(filename):
    tonto.tonto_read_cif(filename.encode("utf-8"), len(filename))

def cleanup():
    tonto.tonto_cleanup()


def get_surface(mol_index):
    n_faces_p = ffi.new("int *")
    n_vertices_p = ffi.new("int *")

    # This actually calculates the HS, after which we can extract
    # the number of vertices/number of faces in the mesh
    tonto.tonto_hs_info(mol_index, n_faces_p, n_vertices_p)
    n_vertices, n_faces = n_vertices_p[0], n_faces_p[0]

    # IMPORTANT TO HAVE ZEROS i.e. allocate the array
    vertices = np.zeros((3, n_vertices), order='F', dtype=np.float64)
    faces = np.zeros((3, n_faces), order='F', dtype=np.intc)

    vertices_p = as_pointer(vertices, "double *")
    faces_p = as_pointer(faces, "int *")

    tonto.tonto_hs_arrays(vertices_p, faces_p, n_vertices, n_faces)
    return vertices, faces

def as_pointer(arr, kind):
    assert arr.flags['F_CONTIGUOUS'], "Array is not fortran contiguous"
    return ffi.cast(kind, arr.__array_interface__['data'][0])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('cif')
    parser.add_argument('-l', '--lmax', default=20, type=int,
                        help='Maximum angular momentum')
    parser.add_argument('-r', '--resolution', default=0.2, type=float,
                        help='Resolution of HS (0.2 = High, 0.5 = Low)')
    parser.add_argument('-b', '--basis', default='/Users/prs/basis_sets',
                        help='Basis set directory to pass to HS program')
    parser.add_argument('--hdf', default='output.h5',
                        help='filename for hdf5 dump file')

    args = parser.parse_args()
    initialize(args)
    read_cif(args.cif)
    for molecule in tqdm(range(1, tonto.n_molecules + 1)):
        vertices, faces = get_surface(molecule)
    tonto.tonto_close_cif();
    read_cif('/Users/prs/repos/crystalexplorer/examples/simple/ETANOL.cif')

    for molecule in tqdm(range(1, tonto.n_molecules + 1)):
        vertices, faces = get_surface(molecule)
    cleanup()

if __name__ == '__main__':
    main()
