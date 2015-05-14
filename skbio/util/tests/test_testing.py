# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import os

import pandas as pd
import numpy.testing as npt

from skbio import OrdinationResults
from skbio.util import get_data_path, assert_ordination_results_equal


def test_get_data_path():
    fn = 'parrot'
    path = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(path, 'data', fn)
    data_path_2 = get_data_path(fn)
    npt.assert_string_equal(data_path_2, data_path)


def test_assert_ordination_results_equal():
    minimal1 = OrdinationResults('foo', 'bar', pd.Series([1.0, 2.0]),
                                 pd.DataFrame([[1, 2, 3], [4, 5, 6]]))

    # a minimal set of results should be equal to itself
    assert_ordination_results_equal(minimal1, minimal1)

    # type mismatch
    with npt.assert_raises(AssertionError):
        assert_ordination_results_equal(minimal1, 'foo')

    # numeric values should be checked that they're almost equal
    almost_minimal1 = OrdinationResults(
        'foo', 'bar',
        pd.Series([1.0000001, 1.9999999]),
        pd.DataFrame([[1, 2, 3], [4, 5, 6]]))
    assert_ordination_results_equal(minimal1, almost_minimal1)

    # test each of the optional numeric attributes
    for attr in ('features', 'samples', 'biplot_scores', 'sample_constraints'):
        # missing optional numeric attribute in one, present in the other
        setattr(almost_minimal1, attr, pd.DataFrame([[1, 2], [3, 4]]))
        with npt.assert_raises(AssertionError):
            assert_ordination_results_equal(minimal1, almost_minimal1)
        setattr(almost_minimal1, attr, None)

        # optional numeric attributes present in both, but not almost equal
        setattr(minimal1, attr, pd.DataFrame([[1, 2], [3, 4]]))
        setattr(almost_minimal1, attr, pd.DataFrame([[1, 2], [3.00002, 4]]))
        with npt.assert_raises(AssertionError):
            assert_ordination_results_equal(minimal1, almost_minimal1)
        setattr(minimal1, attr, None)
        setattr(almost_minimal1, attr, None)

        # optional numeric attributes present in both, and almost equal
        setattr(minimal1, attr, pd.DataFrame([[1.0, 2.0], [3.0, 4.0]]))
        setattr(almost_minimal1, attr,
                pd.DataFrame([[1.0, 2.0], [3.00000002, 4]]))
        assert_ordination_results_equal(minimal1, almost_minimal1)
        setattr(minimal1, attr, None)
        setattr(almost_minimal1, attr, None)

    # missing optional numeric attribute in one, present in the other
    almost_minimal1.proportion_explained = pd.Series([1, 2, 3])
    with npt.assert_raises(AssertionError):
        assert_ordination_results_equal(minimal1, almost_minimal1)
    almost_minimal1.proportion_explained = None

    # optional numeric attributes present in both, but not almost equal
    minimal1.proportion_explained = pd.Series([1, 2, 3])
    almost_minimal1.proportion_explained = pd.Series([1, 2, 3.00002])
    with npt.assert_raises(AssertionError):
        assert_ordination_results_equal(minimal1, almost_minimal1)
    almost_minimal1.proportion_explained = None
    almost_minimal1.proportion_explained = None

    # optional numeric attributes present in both, and almost equal
    minimal1.proportion_explained = pd.Series([1, 2, 3])
    almost_minimal1.proportion_explained = pd.Series([1, 2, 3.00000002])
    with npt.assert_raises(AssertionError):
        assert_ordination_results_equal(minimal1, almost_minimal1)
    almost_minimal1.proportion_explained = None
    almost_minimal1.proportion_explained = None


if __name__ == '__main__':
    import nose
    nose.runmodule()
