# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
from .subsample import subsample, uneven_subsample

__all__ = ['subsample', 'uneven_subsample']

from numpy.testing import Tester
test = Tester().test
