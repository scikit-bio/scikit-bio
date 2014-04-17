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
from nose.tools import assert_almost_equal

from skbio.maths.diversity.alpha.ace import ace


def test_ace():
    assert_almost_equal(ace(np.array([2, 0])), 1.0)
    assert_almost_equal(ace(np.array([12, 0, 9])), 2.0)
    assert_almost_equal(ace(np.array([12, 2, 8])), 3.0)
    assert_almost_equal(ace(np.array([12, 2, 1])), 4.0)
    assert_almost_equal(ace(np.array([12, 1, 2, 1])), 7.0)
    assert_almost_equal(ace(np.array([12, 3, 2, 1])), 4.6)
    assert_almost_equal(ace(np.array([12, 3, 6, 1, 10])), 5.62749672)

    # Just returns the number of species when all are abundant.
    assert_almost_equal(ace(np.array([12, 12, 13, 14])), 4.0)


if __name__ == '__main__':
    import nose
    nose.runmodule()
