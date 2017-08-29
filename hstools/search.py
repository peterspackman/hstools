"""Search module"""
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import os
from pathlib import Path
from scipy.spatial import cKDTree as KDTree
import numpy as np
from collections import namedtuple
import sbf
from .decompose import describe_surface, combination_desc

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

    @staticmethod
    def from_csd_data(l_max=20, use_radius=True):
        """Construct a CSD matcher based on the bundled data

        Keyword arguments:
        l_max -- maximum angular momenta to use for invariants
        (default 20)
        use_radius -- use the mean radius as the first invariant
        (default True)
        """
        names, invariants, radii = load_default_data()
        if use_radius:
            invariants = np.insert(invariants, 0, radii, axis=1)
        return ShapeMatcher(names, invariants)

    @staticmethod
    def from_shapes(shapes, l_max=20):
        """Construct a shapematcher object from a list of shapes

        Arguments:
        shapes -- A list of Shape objects
        Keyword arguments:
        l_max -- maximuma angular momenta to use for invariants
        (default 20)
        """
        invariants, names = [], []
        for name, s in shapes.items():
            invariants.append(s.invariants)
            names.append(name)
        invariants = np.array(invariants)
        names = np.array(names, dtype='|S64')
        return ShapeMatcher(names, invariants)


def chemical_formula(search_result):
    """Convenience function to extract molecular formula from the names/ids
    in the 'universe' dataset

    Arguments:
    search_result -- SearchResult object
    """
    return search.result.name.split('-')[1].split('_')[0]


def surface_description(sbf_file, property_name='shape', properties=None, use_radius=True):
    """Describe a shape/isosurface using spherical harmonics.
    Returns a Shape object

    Arguments:
    sbf_file -- filename or Path object locating a valid surface file

    Keyword Arguments:
    use_radius -- include mean_radius as the first invariant
    """
    log.debug('Describing surface with spherical harmonics')
    if not properties:
        properties = []

    if property_name != 'shape' and property_name not in properties:
        properties.append(property_name)
    if property_name == 'combined':
        name, r, esp, coeffs = combination_desc(sbf_file)
        invariants = make_invariants(coeffs)
        invariants = np.insert(invariants, 0, [r, esp])
    else:
        name, r, coeffs = describe_surface(sbf_file, properties=properties)
        invariants  = make_invariants(coeffs[property_name])
        log.debug('Making invariants')
        if property_name == 'shape' and use_radius:
            invariants = np.insert(invariants, 0, r)
    return Shape(name, invariants)


def load_default_data(directory=os.path.dirname(__file__)):
    """Load the data included with this module.

    Keyword arguments:
    directory -- the folder containing data to load
    (default is the location of this file)
    """
    contents = sbf.read_file(os.path.join(directory, 'main_group.sbf'))
    names = contents['names'].data
    dims = names.shape
    names = names.view('S{}'.format(dims[1])).reshape((dims[0],))
    invariants = contents['invariants'].data
    radii = contents['radii'].data
    return names, invariants, radii


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
        invariants[i] = np.sum(coefficients[l:u+1] *
                               np.conj(coefficients[l:u+1])).real
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
    parser.add_argument('--results', '-n', default=10, type=int,
                        help='Number of matches to print out per surface')
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
    matcher = ShapeMatcher.from_csd_data()
    log.debug('Done')
    paths = list(Path(args.directory).glob(
                '*{}'.format(args.suffix)))
    log.info('%d paths to describe with sht', len(paths))
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(surface_description, path)
                   for path in paths]
        for f in as_completed(futures):
            shape = f.result()
            results = matcher.search_shape(shape, n=args.results)
            print('{}:'.format(shape.name))
            for result in results:
                print('\t{}'.format(result.name))
    log.info('Finished in "%s"', args.directory)

if __name__ == '__main__':
    main()
