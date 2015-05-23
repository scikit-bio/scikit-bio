# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.utils import PY3

import os
import inspect

from nose import core
from nose.tools import nottest

import numpy as np
import numpy.testing as npt
from pandas.util.testing import assert_index_equal


@nottest
class TestRunner(object):
    """Simple wrapper class around nosetests functionality.

    Parameters
    ----------
    filename : str
        __file__ attribute passed in from the caller. This tells the
        tester where to start looking for tests.

    Notes
    -----
    The primary purpose of this class is to create an interface which users
    of scikit-bio can use to run all of the built in tests. Normally this
    would be done by invoking nosetests directly from the command line, but
    scikit-bio needs several additional options which make the command long
    and ugly. This class invokes nose with the required options.

    """
    def __init__(self, filename):
        self._filename = filename
        self._test_dir = os.path.dirname(filename)

    def test(self, verbose=False):
        """Performs the actual running of the tests.

        Parameters
        ----------
        verbose : bool
            flag for running in verbose mode.

        Returns
        -------
        bool
            test run success status
        """
        # NOTE: it doesn't seem to matter what the first element of the argv
        # list is, there just needs to be something there.
        argv = [self._filename, '-I DO_NOT_IGNORE_ANYTHING']
        if not PY3:
            argv.append('--with-doctest')
        if verbose:
            argv.append('-v')
        return core.run(argv=argv, defaultTest=self._test_dir)


def get_data_path(fn, subfolder='data'):
    """Return path to filename ``fn`` in the data folder.

    During testing it is often necessary to load data files. This
    function returns the full path to files in the ``data`` subfolder
    by default.

    Parameters
    ----------
    fn : str
        File name.

    subfolder : str, defaults to ``data``
        Name of the subfolder that contains the data.


    Returns
    -------
    str
        Inferred absolute path to the test data for the module where
        ``get_data_path(fn)`` is called.

    Notes
    -----
    The requested path may not point to an existing file, as its
    existence is not checked.

    """
    # getouterframes returns a list of tuples: the second tuple
    # contains info about the caller, and the second element is its
    # filename
    callers_filename = inspect.getouterframes(inspect.currentframe())[1][1]
    path = os.path.dirname(os.path.abspath(callers_filename))
    data_path = os.path.join(path, subfolder, fn)
    return data_path


def assert_ordination_results_equal(left, right, ignore_method_names=False,
                                    ignore_axis_labels=False,
                                    ignore_biplot_scores_labels=False,
                                    ignore_directionality=False,
                                    decimal=7):
    """Assert that ordination results objects are equal.

    This is a helper function intended to be used in unit tests that need to
    compare ``OrdinationResults`` objects.

    Parameters
    ----------
    left, right : OrdinationResults
        Ordination results to be compared for equality.
    ignore_method_names : bool, optional
        Ignore differences in `short_method_name` and `long_method_name`.
    ignore_axis_labels : bool, optional
        Ignore differences in axis labels (i.e., column labels).
    ignore_biplot_scores_labels : bool, optional
        Ignore differences in `biplot_scores` row and column labels.
    ignore_directionality : bool, optional
        Ignore differences in directionality (i.e., differences in signs) for
        attributes `samples` and `features`.

    Raises
    ------
    AssertionError
        If the two objects are not equal.

    """
    npt.assert_equal(type(left) is type(right), True)

    if not ignore_method_names:
        npt.assert_equal(left.short_method_name, right.short_method_name)
        npt.assert_equal(left.long_method_name, right.long_method_name)

    _assert_frame_equal(left.samples, right.samples,
                        ignore_columns=ignore_axis_labels,
                        ignore_directionality=ignore_directionality,
                        decimal=decimal)

    _assert_frame_equal(left.features, right.features,
                        ignore_columns=ignore_axis_labels,
                        ignore_directionality=ignore_directionality,
                        decimal=decimal)

    _assert_frame_equal(left.biplot_scores, right.biplot_scores,
                        ignore_biplot_scores_labels,
                        ignore_biplot_scores_labels,
                        decimal=decimal)

    _assert_frame_equal(left.sample_constraints, right.sample_constraints,
                        ignore_columns=ignore_axis_labels,
                        decimal=decimal)

    _assert_series_equal(left.eigvals, right.eigvals, ignore_axis_labels,
                         decimal=decimal)

    _assert_series_equal(left.proportion_explained, right.proportion_explained,
                         ignore_axis_labels,
                         decimal=decimal)


def _assert_series_equal(left_s, right_s, ignore_index=False, decimal=7):
    # assert_series_equal doesn't like None...
    if left_s is None or right_s is None:
        assert left_s is None and right_s is None
    else:
        npt.assert_almost_equal(left_s.values, right_s.values,
                                decimal=decimal)
        if not ignore_index:
            assert_index_equal(left_s.index, right_s.index)


def _assert_frame_equal(left_df, right_df, ignore_index=False,
                        ignore_columns=False, ignore_directionality=False,
                        decimal=7):
    # assert_frame_equal doesn't like None...
    if left_df is None or right_df is None:
        assert left_df is None and right_df is None
    else:
        left_values = left_df.values
        right_values = right_df.values

        if ignore_directionality:
            left_values, right_values = _normalize_signs(left_values,
                                                         right_values)
        npt.assert_almost_equal(left_values, right_values, decimal=decimal)

        if not ignore_index:
            assert_index_equal(left_df.index, right_df.index)
        if not ignore_columns:
            assert_index_equal(left_df.columns, right_df.columns)


def _normalize_signs(arr1, arr2):
    """Change column signs so that "column" and "-column" compare equal.

    This is needed because results of eigenproblmes can have signs
    flipped, but they're still right.

    Notes
    =====

    This function tries hard to make sure that, if you find "column"
    and "-column" almost equal, calling a function like np.allclose to
    compare them after calling `normalize_signs` succeeds.

    To do so, it distinguishes two cases for every column:

    - It can be all almost equal to 0 (this includes a column of
      zeros).
    - Otherwise, it has a value that isn't close to 0.

    In the first case, no sign needs to be flipped. I.e., for
    |epsilon| small, np.allclose(-epsilon, 0) is true if and only if
    np.allclose(epsilon, 0) is.

    In the second case, the function finds the number in the column
    whose absolute value is largest. Then, it compares its sign with
    the number found in the same index, but in the other array, and
    flips the sign of the column as needed.
    """
    # Let's convert everyting to floating point numbers (it's
    # reasonable to assume that eigenvectors will already be floating
    # point numbers). This is necessary because np.array(1) /
    # np.array(0) != np.array(1.) / np.array(0.)
    arr1 = np.asarray(arr1, dtype=np.float64)
    arr2 = np.asarray(arr2, dtype=np.float64)

    if arr1.shape != arr2.shape:
        raise ValueError(
            "Arrays must have the same shape ({0} vs {1}).".format(arr1.shape,
                                                                   arr2.shape)
            )

    # To avoid issues around zero, we'll compare signs of the values
    # with highest absolute value
    max_idx = np.abs(arr1).argmax(axis=0)
    max_arr1 = arr1[max_idx, range(arr1.shape[1])]
    max_arr2 = arr2[max_idx, range(arr2.shape[1])]

    sign_arr1 = np.sign(max_arr1)
    sign_arr2 = np.sign(max_arr2)

    # Store current warnings, and ignore division by zero (like 1. /
    # 0.) and invalid operations (like 0. / 0.)
    wrn = np.seterr(invalid='ignore', divide='ignore')
    differences = sign_arr1 / sign_arr2
    # The values in `differences` can be:
    #    1 -> equal signs
    #   -1 -> diff signs
    #   Or nan (0/0), inf (nonzero/0), 0 (0/nonzero)
    np.seterr(**wrn)

    # Now let's deal with cases where `differences != \pm 1`
    special_cases = (~np.isfinite(differences)) | (differences == 0)
    # In any of these cases, the sign of the column doesn't matter, so
    # let's just keep it
    differences[special_cases] = 1

    return arr1 * differences, arr2
