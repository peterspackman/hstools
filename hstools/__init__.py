all = ['decompose', 'lebedev', 'search', 'gui']
from .search import get_csd_matcher
from .search import ShapeMatcher, Shape
from .search import get_shape_description as describe
from warnings import warn
import numpy as np
__MATCHER__ = None

class MatcherBox():
    pass

__matcher = MatcherBox()
__matcher.universe_matcher = None

def csd_matcher():
    if __matcher.universe_matcher is None:
        warn('Initializing matcher, this make take a few seconds')
        __matcher.universe_matcher = get_csd_matcher()
    return __matcher.universe_matcher
