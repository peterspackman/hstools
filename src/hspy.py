#! /usr/bin/env python3
import os
import sys
from cffi import FFI
import numpy as np
import argparse

ffi = FFI()
ffi.cdef('''
void f_print_state();
void f_set_basis_directory(const char *input_string, const int length);
void f_set_cif(const char *input_string, const int length);
void f_set_lmax(int l);
void f_set_resolution(double r);
void f_initialize();
void f_cleanup();
void f_make_surfaces();
''')

if sys.platform == 'darwin':
    suffix = 'dylib'
else:
    suffx = 'so'

tonto = ffi.dlopen('src/libhsapi.{}'.format(suffix))

def initialize(args):
    tonto.f_set_basis_directory(args.basis.encode("utf-8"),
                                len(args.basis))
    tonto.f_set_cif(args.cif.encode("utf-8"), len(args.cif))
    tonto.f_set_lmax(args.lmax)
    tonto.f_set_resolution(args.resolution)
    tonto.f_initialize()

def cleanup():
    tonto.f_cleanup()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('cif')
    parser.add_argument('-l', '--lmax', default=20, type=int,
                        help='Maximum angular momentum')
    parser.add_argument('-r', '--resolution', default=0.2, type=float,
                        help='Resolution of HS (0.2 = High, 0.5 = Low)')
    parser.add_argument('-b', '--basis', default='/Users/prs/basis_sets',
                        help='Basis set directory to pass to HS program')
    args = parser.parse_args()
    initialize(args)
    tonto.f_make_surfaces()
    cleanup()

if __name__ == '__main__':
    main()
