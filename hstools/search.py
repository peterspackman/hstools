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
        """Match other shapes based on euclidean distance. 
        Constructs a KDTree in order to do nearest neighbour queries.
        For large datasets it might take a second or two to build the tree.
        
        Arguments:
        ids -- set names/identifiers for the shapes
        invariants -- 2D array of invariants that describe the shapes
        """
        self.ids = ids
        self.invariants = invariants
        log.debug('Constructing tree from %d invariants', len(invariants))
        self.tree = KDTree(invariants)

    def search_invariants(self, invariants, n=10):
        """Search for matches based on invariants.

        Arguments:
        invariants -- N length array of shape descriptors

        Keyword arguments:
        n -- number of matches to return (default 10)
        """
        log.debug('Searching for %d closest points', n)
        distances, indexes = self.tree.query(invariants, n)
        return [SearchResult(name.decode('utf-8'), d, inv) 
                for name, d, inv in 
                zip(self.ids[indexes], distances, self.invariants[indexes])]

    def search_shape(self, shape, n=10):
        """Search for matches based on a shape object. (convenience function)

        Arguments:
        shape -- a Shape object.

        Keyword arguments:
        n -- number of matches to return (default 10)
        """
        log.debug('Searching for closest shapes to %s', shape.name)
        return self.search_invariants(shape.invariants, n=n)


def get_chemical_formula(search_result):
    """Convenience function to extract molecular formula from the names/ids
    in the 'universe' dataset

    Arguments:
    search_result -- SearchResult object
    """
    return search.result.name.split('-')[1].split('_')[0]


def get_shape_description(sbf_file, use_radius=True):
    """Describe a shape/isosurface using spherical harmonics.
    Returns a Shape object

    Arguments:
    sbf_file -- filename or Path object locating a valid surface file

    Keyword Arguments:
    use_radius -- include mean_radius as the first invariant
    """
    log.debug('Describing surface with spherical harmonics')
    name, r, coeffs = describe_surface(sbf_file)
    log.debug('Making invariants')
    invariants = np.insert(make_invariants(coeffs), 0, r)
    return Shape(name, invariants)
    

def load_default_data(directory=os.path.dirname(__file__)): 
    """Load the data included with this module.
    
    Keyword arguments:
    directory -- the folder containing data to load 
    (default is the location of this file)
    """
    names = np.fromfile(os.path.join(directory, 'universe-names.bin'), dtype='S64')
    invariants = np.fromfile(
            os.path.join(directory, 'universe-invariants.bin'),
            dtype=np.float64).reshape(names.size, 21)

    radii = np.fromfile(
            os.path.join(directory, 'universe-radii.bin'),
            dtype=np.float64)
    
    return names, invariants, radii


def get_csd_matcher(l_max=20, use_radius=True):
    names, invariants, radii = load_default_data()
    if use_radius:
        invariants = np.insert(invariants, 0, radii, axis=1)
    return ShapeMatcher(names, invariants)


def make_invariants(coefficients):
    """Construct the 'N' type invariants from sht coefficients.
    If coefficients is of length n, the size of the result will be sqrt(n)

    Arguments:
    coefficients -- the set of spherical harmonic coefficients
    """
    size = int(np.sqrt(len(coefficients)))
    invariants = np.empty(shape=(size), dtype=np.float64)
    for i in range(0, size):
        l, u = i**2, (i+1)**2
        invariants[i] = np.sum(coefficients[l:u+1] * np.conj(coefficients[l:u+1])).real
    return invariants

def add_files_from_directory(directory, data_dict={}):
    for f in Path(directory).glob('*.npy'):
        data_dict[f.stem] = np.load(f)
    return data_dict


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
    """Read through all sbf files in a directory, writing
    the names of the top matches from the database
    """
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument('directory')
    parser.add_argument('--log-file', default=None,
                        help='Log to file instead of stdout')
    parser.add_argument('--suffix', '-s', default='.sbf',
                        help='File suffix to find sbf files')
    parser.add_argument('--log-level', default='INFO',
                        help='Log level')
    parser.add_argument('--jobs', '-j', default=4, type=int,
                        help='Number of parallel jobs to run')
    args = parser.parse_args()
    if args.log_file:
        logging.basicConfig(filename=args.log_file, level=args.log_level)
    else:
        logging.basicConfig(level=args.log_level)
    log.info('Starting %s, jobs: %d', args.directory, args.jobs)
    log.debug('Creating CSD matcher')
    matcher = get_csd_matcher()
    log.debug('Done')
    paths = list(Path(args.directory).glob('*{}'.format(args.suffix)))
    log.info('%d paths to describe with sht', len(paths))
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(get_shape_description, path) for path in paths]
        for f in as_completed(futures):
            shape = f.result()
            results = matcher.search_shape(shape)
            print('{}:'.format(shape.name))
            for result in results:
                print('\t{}'.format(result.name))
    log.info('Finished in "%s"', args.directory)

if __name__ == '__main__':
    main()
