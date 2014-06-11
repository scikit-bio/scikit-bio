#!/usr/bin/env python
from __future__ import absolute_import, division, print_function
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os

import numpy.testing as npt

from skbio.util.testing import get_data_path


def test_get_data_path():
    fn = 'parrot'
    path = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(path, 'data', fn)
    data_path_2 = get_data_path(fn)
    npt.assert_string_equal(data_path_2, data_path)


if __name__ == '__main__':
    import nose
    nose.runmodule()
