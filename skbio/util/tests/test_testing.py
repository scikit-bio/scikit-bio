# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import os
import itertools
import unittest

import pandas as pd
import numpy as np

from skbio.util import get_data_path, assert_data_frame_almost_equal


class TestGetDataPath(unittest.TestCase):
    def test_get_data_path(self):
        fn = 'parrot'
        path = os.path.dirname(os.path.abspath(__file__))
        data_path = os.path.join(path, 'data', fn)
        data_path_2 = get_data_path(fn)
        self.assertEqual(data_path_2, data_path)


class TestAssertDataFrameAlmostEqual(unittest.TestCase):
    def setUp(self):
        self.df = pd.DataFrame(
            {'foo': [42, 42.0, np.nan, 0],
             'bar': ['a', 'b', 'cd', 'e']})

    def test_not_equal(self):
        unequal_dfs = [
            self.df,
            # floating point error too large to be "almost equal"
            pd.DataFrame({'foo': [42, 42.001, np.nan, 0],
                          'bar': ['a', 'b', 'cd', 'e']}),
            # extra NaN
            pd.DataFrame({'foo': [42, np.nan, np.nan, 0],
                          'bar': ['a', 'b', 'cd', 'e']}),
            # different column order
            pd.DataFrame(self.df, columns=['foo', 'bar']),
            # different index order
            pd.DataFrame(self.df, index=np.arange(4)[::-1]),
            # different index type
            pd.DataFrame(self.df, index=np.arange(4).astype(float)),
            # various forms of "empty" DataFrames that are not equivalent
            pd.DataFrame(),
            pd.DataFrame(index=np.arange(10)),
            pd.DataFrame(columns=np.arange(10)),
            pd.DataFrame(index=np.arange(10), columns=np.arange(10)),
            pd.DataFrame(index=np.arange(9)),
            pd.DataFrame(columns=np.arange(9)),
            pd.DataFrame(index=np.arange(9), columns=np.arange(9))
        ]

        # each df should compare equal to itself
        for df in unequal_dfs:
            assert_data_frame_almost_equal(df, df)

        # every pair of dfs should not compare equal. use permutations instead
        # of combinations to test that comparing df1 to df2 and df2 to df1 are
        # both not equal
        for df1, df2 in itertools.permutations(unequal_dfs, 2):
            with self.assertRaises(AssertionError):
                assert_data_frame_almost_equal(df1, df2)

    def test_equal(self):
        equal_dfs = [
            self.df,
            # floating point error small enough to be "almost equal"
            pd.DataFrame({'foo': [42, 42.00001, np.nan, 0],
                          'bar': ['a', 'b', 'cd', 'e']})
        ]

        for df in equal_dfs:
            assert_data_frame_almost_equal(df, df)

        for df1, df2 in itertools.permutations(equal_dfs, 2):
            assert_data_frame_almost_equal(df1, df2)


if __name__ == '__main__':
    unittest.main()
