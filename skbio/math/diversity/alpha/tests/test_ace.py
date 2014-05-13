#!/usr/bin/env python
from __future__ import division

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
from nose.tools import assert_almost_equal, assert_raises, assert_true

from skbio.math.diversity.alpha.ace import ace


def test_ace():
    assert_almost_equal(ace(np.array([2, 0])), 1.0)
    assert_almost_equal(ace(np.array([12, 0, 9])), 2.0)
    assert_almost_equal(ace(np.array([12, 2, 8])), 3.0)
    assert_almost_equal(ace(np.array([12, 2, 1])), 4.0)
    assert_almost_equal(ace(np.array([12, 1, 2, 1])), 7.0)
    assert_almost_equal(ace(np.array([12, 3, 2, 1])), 4.6)
    assert_almost_equal(ace(np.array([12, 3, 6, 1, 10])), 5.62749672)

    # Just returns the number of OTUs when all are abundant.
    assert_almost_equal(ace(np.array([12, 12, 13, 14])), 4.0)

    # Border case: only singletons and 10-tons, no abundant OTUs.
    assert_almost_equal(ace([0, 1, 1, 0, 0, 10, 10, 1, 0, 0]), 9.35681818182)


def test_ace_only_rare_singletons():
    with assert_raises(ValueError):
        ace([0, 0, 43, 0, 1, 0, 1, 42, 1, 43])


if __name__ == '__main__':
    import nose
    nose.runmodule()
