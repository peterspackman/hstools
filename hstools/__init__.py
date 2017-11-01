from .search import ShapeMatcher, Shape
import os
all = ['decompose', 'lebedev', 'search', 'gui', 'notebook']
__MATCHER__ = None


class UnknownMatcherException(Exception):
    pass

class MatcherBox():
    """Simple storage class for caching
    data/matcher"""
    pass

__matcher = MatcherBox()
__matcher.matcher_hs = None
__matcher.matcher_ps = None
__data_location = os.path.dirname(__file__)
HS_DATA_LOCATION = os.path.join(__data_location, 'shapes-hs.sbf')
PS_DATA_LOCATION = os.path.join(__data_location, 'shapes-ps.sbf')


def csd_matcher(kind='hirshfeld'):
    """Create the matcher from the bundled data extracted
    from the CSD"""
    if kind == 'hirshfeld':
        if __matcher.matcher_hs is None:
            __matcher.matcher_hs = ShapeMatcher.from_datafile(HS_DATA_LOCATION)
        return __matcher.matcher_hs
    elif kind == 'promolecule':
        if __matcher.matcher_ps is None:
            __matcher.matcher_ps = ShapeMatcher.from_datafile(PS_DATA_LOCATION)
        return __matcher.matcher_ps
    else:
        raise UnknownMatcherException(kind)

