from scipy.spatial import KDTree
import numpy as np
import logging
from pathlib import Path
from collections import namedtuple
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from .decompose import describe_surface

log = logging.getLogger(__name__)
Shape = namedtuple('Shape', 'name invariants')
SearchResult = namedtuple('SearchResult', 'name proximity invariants')

class ShapeMatcher(object):
    def __init__(self, ids, invariants):
        self.ids = ids
        self.invariants = invariants
        log.debug('Constructing tree from %d invariants', len(invariants))
        self.tree = KDTree(invariants)

    def search_invariants(self, invariants, n=10):
        log.debug('Searching for %d closest points', n)
        distances, indexes = self.tree.query(invariants, n)
        return [SearchResult(name.decode('utf-8'), d, inv) 
                for name, d, inv in 
                zip(self.ids[indexes], distances, self.invariants[indexes])]

    def search_shape(self, shape, n=10):
        log.debug('Searching for closest shapes to %s', shape.name)
        return self.search_invariants(shape.invariants, n=n)


def get_chemical_formula(search_result):
    return search.result.name.split('-')[1].split('_')[0]


def get_shape_description(sbf_file, use_radius=True):
    log.debug('Describing surface with spherical harmonics')
    name, r, coeffs = describe_surface(sbf_file)
    log.debug('Making invariants')
    invariants = np.insert(make_invariants(coeffs), 0, r)
    return Shape(name, invariants)
    

def load_default_data(directory=os.path.dirname(__file__)): 
    names = np.fromfile(os.path.join(directory, 'universe-names.bin'), dtype='S64')
    #coefficients = np.fromfile(os.path.join(directory, 'universe-coeffs.bin'))
    invariants = np.fromfile(
            os.path.join(directory, 'universe-invariants.bin'),
            dtype=np.float64).reshape(names.size, 21)

    radii = np.fromfile(
            os.path.join(directory, 'universe-radii.bin'),
            dtype=np.float64)
    
    return names, invariants, radii

def add_files_from_directory(directory, data_dict={}):
    for f in Path(directory).glob('*.npy'):
        data_dict[f.stem] = np.load(f)
    return data_dict

def make_invariants(coefficients):
    size = 20 + 1
    invariants = np.empty(shape=(size), dtype=np.float64)
    for i in range(0, size):
        l, u = i**2, (i+1)**2
        invariants[i] = np.sum(coefficients[l:u+1] * np.conj(coefficients[l:u+1])).real
    return invariants

def create_arrays(data_dict):
    n = len(data_dict)
    names = np.empty(shape=(n), dtype='S64')
    coefficients = np.empty(shape=(n, 441), dtype=np.complex128)
    invariants = np.empty(shape=(n, 20+1), dtype=np.float64)
    keys, values = data_dict.keys(), data_dict.values()
    names[:] = [key for key in keys]
    coefficients[:, :] = [value for value in values]
    invariants[:] = [make_invariants(c) for c in values]
    return names, coefficients, invariants

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
    
