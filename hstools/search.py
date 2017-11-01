"""Search module"""
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import os
from pathlib import Path
from scipy.spatial import cKDTree as KDTree
import numpy as np
from collections import namedtuple
from .decompose import sht_isosurface, Shape
import pandas as pd
import sbf
import shtns

LOG = logging.getLogger(__name__)
__SearchResult = namedtuple('SearchResult', 'name proximity invariants')

class SearchResult(__SearchResult):
    @property
    def chemical_formula(self):
        """Convenience function to extract molecular formula from the names/ids
        in the provided dataset
        """
        return self.name.split('-')[1].split('_')[0]

    @property
    def csd_refcode(self):
        """Convenience function to extract CSD refcode from the names/ids
        in the provided datasets
        """
        return self.name.split('-')[0]


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
        LOG.debug('Constructing tree from %d invariants', len(invariants))
        self.tree = KDTree(invariants)

    def search_invariants(self, invariants, n=10, df=False):
        """Search for matches based on invariants.

        Arguments:
        invariants -- N length array of shape descriptors

        Keyword arguments:
        n -- number of matches to return (default 10)
        df -- return matches as a pandas DataFrame (default False)
        """
        if n == 'max':
            n = len(self.invariants)
        LOG.debug('Searching for %d closest points', n)
        distances, indexes = self.tree.query(invariants, n)
        invariants = self.invariants[indexes]
        # Need to handle case of n == 1 correctly
        if isinstance(indexes, int):
            ids = self.ids[indexes].decode('utf-8')
            return SearchResult(ids, distances, invariants)
        else:
            ids = [x.decode('utf-8') for x in self.ids[indexes]]
        if df:
            return pd.DataFrame({'ID': ids, 'Proximity': distances}).set_index('ID')
        else:
            return [SearchResult(n, d, i)
                        for n, d, i in
                        zip(ids, distances, invariants)]

    def search_shape(self, shape, **kwargs):
        """Search for matches based on a shape object. (convenience function)

        Arguments:
        shape -- a Shape object.

        Keyword arguments:
        n -- number of matches to return (default 10)
        df -- return matches as a pandas DataFrame (default False)
        """
        LOG.debug('Searching for closest shapes to %s', shape.name)
        # delegate to search_invariants method
        return self.search_invariants(shape.invariants, **kwargs)

    @staticmethod
    def from_datafile(filename, l_max=20):
        """Construct a CSD matcher based on the bundled data

        Keyword arguments:
        l_max -- maximum angular momenta to use for invariants
        (default 20)
        use_radius -- use the mean radius as the first invariant
        (default True)
        """
        names, invariants = load_data(filename)
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
        if isinstance(shapes, dict):
            for name, s in shapes.items():
                invariants.append(s.invariants)
                names.append(name)
        else:
            for s in shapes:
                invariants.append(s.invariants)
                names.append(s.name)
        invariants = np.array(invariants)
        names = np.array(names, dtype='|S64')
        return ShapeMatcher(names, invariants)

    @staticmethod
    def from_surface_files(files, property_name='shape'):
        """Construct a CSD matcher based on the bundled data

        Keyword arguments:
        l_max -- maximum angular momenta to use for invariants
        (default 20)
        use_radius -- use the mean radius as the first invariant
        (default True)
        """
        shapes = {}
        for f in files:
            shapes['f.stem'] = surface_description(f, property_name=property_name)

        return ShapeMatcher.from_shapes(shapes)





def load_data(filename):
    """Load the data included with this module.

    Keyword arguments:
    directory -- the folder containing data to load
    (default is the location of this file)
    """
    contents = sbf.read_file(filename)
    names = contents['names'].data
    dims = names.shape
    names = names.view('S{}'.format(dims[1])).reshape((dims[0],))
    invariants = contents['invariants'].data
    return names, invariants



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
    parser.add_argument('--suffix', '-s', default='-hs.sbf',
                        help='File suffix to find sbf files')
    parser.add_argument('--results', '-n', default=10, type=int,
                        help='Number of matches to print out per surface')
    parser.add_argument('--surface-type', '-t', default='hirshfeld', type=str,
                        choices={'hirshfeld', 'promolecule'},
                        help='Surface type to match against')
    parser.add_argument('--log-level', default='INFO',
                        help='Log level')
    parser.add_argument('--jobs', '-j', default=1, type=int,
                        help='Number of parallel jobs to run')
    args = parser.parse_args()
    if args.log_file:
        logging.basicConfig(filename=args.log_file, level=args.log_level)
    else:
        logging.basicConfig(level=args.log_level)
    LOG.info('Starting %s, jobs: %d', args.directory, args.jobs)
    LOG.debug('Creating CSD matcher')
    from . import csd_matcher
    from .decompose import surface_description
    matcher = csd_matcher(args.surface_type)
    LOG.debug('Done')
    paths = list(Path(args.directory).glob(
                '*{}'.format(args.suffix)))
    LOG.info('%d paths to describe with sht', len(paths))
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(surface_description, path)
                   for path in paths]
        for f in as_completed(futures):
            shape = f.result()
            results = matcher.search_shape(shape, n=args.results, df=True)
            print()
            print(shape.name)
            print()
            print(results)
            print()
            print()
    LOG.info('Finished in "%s"', args.directory)

if __name__ == '__main__':
    main()
