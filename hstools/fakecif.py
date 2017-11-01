#!/usr/bin/env python3
import argparse
from pathlib import Path
from collections import namedtuple, defaultdict
import logging
log = logging.getLogger('fakecif')

Atom = namedtuple('Atom', 'element label center')

cif_contents = """# Fake CIF generated from {name}
data_{name}
_symmetry_cell_setting           triclinic
_symmetry_space_group_name_H-M   'P 1'
_symmetry_Int_Tables_number      1
_space_group_name_Hall           'P 1'
loop_
_symmetry_equiv_pos_site_id
_symmetry_equiv_pos_as_xyz
1 x,y,z
_cell_length_a                   {a}(1)
_cell_length_b                   {b}(1)
_cell_length_c                   {c}(1)
_cell_angle_alpha                90.00(1)
_cell_angle_beta                 90.00(1)
_cell_angle_gamma                90.00(1)
_cell_volume                     {cell_volume}
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
{cell_atoms}
"""


def center_string(coord):
    return '\t'.join('{:.8f}'.format(c) for c in coord)


def atom_string(atom):
    return '\t'.join([atom.label, atom.element, center_string(atom.center)])


def output_fake_cif(path, atoms, cell_dims):
    cell_volume = cell_dims[0] * cell_dims[1] * cell_dims[2]
    cell_atoms = '\n'.join(atom_string(atom) for atom in atoms)
    output = cif_contents.format(name=path.stem,
                                 a=cell_dims[0],
                                 b=cell_dims[1],
                                 c=cell_dims[2],
                                 cell_volume=cell_volume,
                                 cell_atoms=cell_atoms)
    with path.open('w') as f:
        f.write(output)


def bounding_box(atoms):
    bounds = []
    for i in range(0, 3):
        minimum = min(a.center[i] for a in atoms)
        maximum = max(a.center[i] for a in atoms)
        bounds.append([minimum, maximum])
    return bounds


def convert_to_fractional_coords(atoms, cell_dims):
    log.debug('Converting to fraction coords')
    for atom in atoms:
        for i in range(0, 3):
            atom.center[i] /= cell_dims[i]


def process_xyz_file(path):
    element_count = defaultdict(int)
    atoms = []
    with path.open() as f:
        n = int(f.readline().strip())
        comment = f.readline().strip()
        for line in f:
            tokens = line.strip().split()
            if (len(tokens) == 4):
                element = tokens[0]
                center = [float(x) for x in tokens[1:4]]
                element_count[element] += 1
                atoms.append(Atom(element,
                                  '{}{}'.format(element,
                                                element_count[element]),
                                  center))
            else:
                print('too many tokens ({})'.format(len(tokens)))
    bounds = bounding_box(atoms)
    cell_dims = [round(x2 - x1) * 1.5 for x1, x2 in bounds]
    dim = max(cell_dims)
# shift the minimum to be at the origin
    for atom in atoms:
        for i in range(0, 3):
            atom.center[i] -= bounds[i][0]
    cell_dims = [dim, dim, dim]
    convert_to_fractional_coords(atoms, cell_dims)
    return atoms, cell_dims


def make_cif(input_path):
    atoms, cell_dims = process_xyz_file(input_path)
    cif_path = input_path.with_suffix('.cif')
    log.info('Writing CIF to %s', cif_path)
    output_fake_cif(cif_path, atoms, cell_dims)


def main():
    """Read through all xyz files in a directory, writing
    fake CIFs for each of them
    """
    import argparse
    import os
    from concurrent.futures import ThreadPoolExecutor, as_completed
    parser = argparse.ArgumentParser()
    parser.add_argument('directory')
    parser.add_argument('--log-file', default=None,
                        help='Log to file instead of stdout')
    parser.add_argument('--suffix', '-s', default='.xyz',
                        help='File suffix to find xyz files')
    parser.add_argument('--log-level', default='INFO',
                        help='Log level')
    parser.add_argument('--jobs', '-j', default=1, type=int,
                        help='Number of parallel jobs to run')
    args = parser.parse_args()

    if args.log_file:
        logging.basicConfig(filename=args.log_file, level=args.log_level)
    else:
        logging.basicConfig(level=args.log_level)

    log.info('Searching in %s', 
             args.directory)

    paths = list(Path(args.directory).glob('*' + args.suffix))
    log.info('%d xyz files to process', len(paths))
    for path in paths:
        make_cif(path)
    log.info('Finished %s', args.directory)
