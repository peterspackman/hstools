#!/usr/bin/env python3
import argparse
from pathlib import Path
from collections import namedtuple, defaultdict

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
_cell_size_a                   {a}(1)
_cell_size_b                   {b}(1)
_cell_size_c                   {c}(1)
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
                                 cell_volume=cell_volume, cell_atoms=cell_atoms)
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
            tokens = line.split()
            if (len(tokens) == 4):
                element = tokens[0]
                center = [float(x) for x in tokens[1:4]]
                element_count[element] += 1
                atoms.append(Atom(element, '{}{}'.format(element, element_count[element]),center))
    bounds = bounding_box(atoms)
    cell_dims = [round(x2 - x1) * 1.05 for x1, x2 in bounds ]
# shift the minimum to be at the origin
    for atom in atoms:
        for i in range(0, 3):
            atom.center[i] -= bounds[i][0]
    return atoms, cell_dims


def process_xyz_files():
    pattern = '*xyz'
    current_directory = Path('.')
    for path in current_directory.glob(pattern):
        atoms, cell_dims = process_xyz_file(path)
        output_fake_cif(path.with_suffix('.cif'), atoms, cell_dims)


if __name__ == '__main__':
    process_xyz_files()

