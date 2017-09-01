from .search import ShapeMatcher, Shape
from warnings import warn
all = ['decompose', 'lebedev', 'search', 'gui', 'notebook']
__MATCHER__ = None


class MatcherBox():
    """Simple storage class for caching
    data/matcher"""
    pass

__matcher = MatcherBox()
__matcher.universe_matcher = None


def csd_matcher():
    """Create the matcher from the bundled data extracted
    from the CSD"""
    if __matcher.universe_matcher is None:
        __matcher.universe_matcher = ShapeMatcher.from_csd_data()
    return __matcher.universe_matcher
