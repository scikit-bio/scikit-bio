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

import six
import pandas as pd
import pandas.util.testing as pdt
from nose import core
from nose.tools import nottest

import numpy as np
import numpy.testing as npt
from pandas.util.testing import assert_index_equal

from ._decorator import experimental


class ReallyEqualMixin(object):
    """Use this for testing __eq__/__ne__.

    Taken and modified from the following public domain code:
      https://ludios.org/testing-your-eq-ne-cmp/

    """

    def assertReallyEqual(self, a, b):
        # assertEqual first, because it will have a good message if the
        # assertion fails.
        self.assertEqual(a, b)
        self.assertEqual(b, a)
        self.assertTrue(a == b)
        self.assertTrue(b == a)
        self.assertFalse(a != b)
        self.assertFalse(b != a)

        # We do not support cmp/__cmp__ because they do not exist in Python 3.
        # However, we still test this to catch potential bugs where the
        # object's parent class defines a __cmp__.
        if not PY3:
            self.assertEqual(0, cmp(a, b))  # noqa
            self.assertEqual(0, cmp(b, a))  # noqa

    def assertReallyNotEqual(self, a, b):
        # assertNotEqual first, because it will have a good message if the
        # assertion fails.
        self.assertNotEqual(a, b)
        self.assertNotEqual(b, a)
        self.assertFalse(a == b)
        self.assertFalse(b == a)
        self.assertTrue(a != b)
        self.assertTrue(b != a)

        # We do not support cmp/__cmp__ because they do not exist in Python 3.
        # However, we still test this to catch potential bugs where the
        # object's parent class defines a __cmp__.
        if not PY3:
            self.assertNotEqual(0, cmp(a, b))  # noqa
            self.assertNotEqual(0, cmp(b, a))  # noqa


class MetadataMixinTests(object):
    def test_constructor_invalid_type(self):
        for md in (0, 'a', ('f', 'o', 'o'), np.array([]), pd.DataFrame()):
            with six.assertRaisesRegex(self, TypeError,
                                       'metadata must be a dict'):
                self._metadata_constructor_(metadata=md)

    def test_constructor_no_metadata(self):
        for md in None, {}:
            obj = self._metadata_constructor_(metadata=md)

            self.assertFalse(obj.has_metadata())
            self.assertEqual(obj.metadata, {})

    def test_constructor_with_metadata(self):
        obj = self._metadata_constructor_(metadata={'foo': 'bar'})
        self.assertEqual(obj.metadata, {'foo': 'bar'})

        obj = self._metadata_constructor_(
                metadata={'': '', 123: {'a': 'b', 'c': 'd'}})
        self.assertEqual(obj.metadata, {'': '', 123: {'a': 'b', 'c': 'd'}})

    def test_constructor_handles_missing_metadata_efficiently(self):
        self.assertIsNone(self._metadata_constructor_()._metadata)
        self.assertIsNone(self._metadata_constructor_(metadata=None)._metadata)

    def test_constructor_makes_shallow_copy_of_metadata(self):
        md = {'foo': 'bar', 42: []}
        obj = self._metadata_constructor_(metadata=md)

        self.assertEqual(obj.metadata, md)
        self.assertIsNot(obj.metadata, md)

        md['foo'] = 'baz'
        self.assertEqual(obj.metadata, {'foo': 'bar', 42: []})

        md[42].append(True)
        self.assertEqual(obj.metadata, {'foo': 'bar', 42: [True]})

    def test_eq(self):
        self.assertReallyEqual(
                self._metadata_constructor_(metadata={'foo': 42}),
                self._metadata_constructor_(metadata={'foo': 42}))

        self.assertReallyEqual(
                self._metadata_constructor_(metadata={'foo': 42, 123: {}}),
                self._metadata_constructor_(metadata={'foo': 42, 123: {}}))

    def test_eq_missing_metadata(self):
        self.assertReallyEqual(self._metadata_constructor_(),
                               self._metadata_constructor_())
        self.assertReallyEqual(self._metadata_constructor_(),
                               self._metadata_constructor_(metadata={}))
        self.assertReallyEqual(self._metadata_constructor_(metadata={}),
                               self._metadata_constructor_(metadata={}))

    def test_eq_handles_missing_metadata_efficiently(self):
        obj1 = self._metadata_constructor_()
        obj2 = self._metadata_constructor_()
        self.assertReallyEqual(obj1, obj2)

        self.assertIsNone(obj1._metadata)
        self.assertIsNone(obj2._metadata)

    def test_ne(self):
        # Both have metadata.
        obj1 = self._metadata_constructor_(metadata={'id': 'foo'})
        obj2 = self._metadata_constructor_(metadata={'id': 'bar'})
        self.assertReallyNotEqual(obj1, obj2)

        # One has metadata.
        obj1 = self._metadata_constructor_(metadata={'id': 'foo'})
        obj2 = self._metadata_constructor_()
        self.assertReallyNotEqual(obj1, obj2)

    def test_metadata_getter(self):
        obj = self._metadata_constructor_(
                metadata={42: 'foo', ('hello', 'world'): 43})

        self.assertIsInstance(obj.metadata, dict)
        self.assertEqual(obj.metadata, {42: 'foo', ('hello', 'world'): 43})

        obj.metadata[42] = 'bar'
        self.assertEqual(obj.metadata, {42: 'bar', ('hello', 'world'): 43})

    def test_metadata_getter_no_metadata(self):
        obj = self._metadata_constructor_()

        self.assertIsNone(obj._metadata)
        self.assertIsInstance(obj.metadata, dict)
        self.assertEqual(obj.metadata, {})
        self.assertIsNotNone(obj._metadata)

    def test_metadata_setter(self):
        obj = self._metadata_constructor_()

        self.assertFalse(obj.has_metadata())

        obj.metadata = {'hello': 'world'}
        self.assertTrue(obj.has_metadata())
        self.assertEqual(obj.metadata, {'hello': 'world'})

        obj.metadata = {}
        self.assertFalse(obj.has_metadata())
        self.assertEqual(obj.metadata, {})

    def test_metadata_setter_makes_shallow_copy(self):
        obj = self._metadata_constructor_()

        md = {'foo': 'bar', 42: []}
        obj.metadata = md

        self.assertEqual(obj.metadata, md)
        self.assertIsNot(obj.metadata, md)

        md['foo'] = 'baz'
        self.assertEqual(obj.metadata, {'foo': 'bar', 42: []})

        md[42].append(True)
        self.assertEqual(obj.metadata, {'foo': 'bar', 42: [True]})

    def test_metadata_setter_invalid_type(self):
        obj = self._metadata_constructor_(metadata={123: 456})

        for md in (None, 0, 'a', ('f', 'o', 'o'), np.array([]),
                   pd.DataFrame()):
            with six.assertRaisesRegex(self, TypeError,
                                       'metadata must be a dict'):
                obj.metadata = md
            self.assertEqual(obj.metadata, {123: 456})

    def test_metadata_deleter(self):
        obj = self._metadata_constructor_(metadata={'foo': 'bar'})

        self.assertEqual(obj.metadata, {'foo': 'bar'})

        del obj.metadata
        self.assertIsNone(obj._metadata)
        self.assertFalse(obj.has_metadata())

        # Delete again.
        del obj.metadata
        self.assertIsNone(obj._metadata)
        self.assertFalse(obj.has_metadata())

        obj = self._metadata_constructor_()

        self.assertIsNone(obj._metadata)
        self.assertFalse(obj.has_metadata())
        del obj.metadata
        self.assertIsNone(obj._metadata)
        self.assertFalse(obj.has_metadata())

    def test_has_metadata(self):
        obj = self._metadata_constructor_()

        self.assertFalse(obj.has_metadata())
        # Handles metadata efficiently.
        self.assertIsNone(obj._metadata)

        self.assertFalse(
                self._metadata_constructor_(metadata={}).has_metadata())

        self.assertTrue(
                self._metadata_constructor_(metadata={'': ''}).has_metadata())
        self.assertTrue(
                self._metadata_constructor_(
                        metadata={'foo': 42}).has_metadata())


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
    @experimental(as_of="0.4.0")
    def __init__(self, filename):
        self._filename = filename
        self._test_dir = os.path.dirname(filename)

    @experimental(as_of="0.4.0")
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
        if PY3:
            argv.extend(['--with-doctest', '--doctest-tests'])
        if verbose:
            argv.append('-v')
        return core.run(argv=argv, defaultTest=self._test_dir)


@experimental(as_of="0.4.0")
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


@experimental(as_of="0.4.0")
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


@experimental(as_of="0.4.0")
def assert_data_frame_almost_equal(left, right):
    """Raise AssertionError if ``pd.DataFrame`` objects are not "almost equal".

    Wrapper of ``pd.util.testing.assert_frame_equal``. Floating point values
    are considered "almost equal" if they are within a threshold defined by
    ``assert_frame_equal``. This wrapper uses a number of
    checks that are turned off by default in ``assert_frame_equal`` in order to
    perform stricter comparisons (for example, ensuring the index and column
    types are the same). It also does not consider empty ``pd.DataFrame``
    objects equal if they have a different index.

    Other notes:

    * Index (row) and column ordering must be the same for objects to be equal.
    * NaNs (``np.nan``) in the same locations are considered equal.

    This is a helper function intended to be used in unit tests that need to
    compare ``pd.DataFrame`` objects.

    Parameters
    ----------
    left, right : pd.DataFrame
        ``pd.DataFrame`` objects to compare.

    Raises
    ------
    AssertionError
        If `left` and `right` are not "almost equal".

    See Also
    --------
    pandas.util.testing.assert_frame_equal

    """
    # pass all kwargs to ensure this function has consistent behavior even if
    # `assert_frame_equal`'s defaults change
    pdt.assert_frame_equal(left, right,
                           check_dtype=True,
                           check_index_type=True,
                           check_column_type=True,
                           check_frame_type=True,
                           check_less_precise=False,
                           check_names=True,
                           by_blocks=False,
                           check_exact=False)
    # this check ensures that empty DataFrames with different indices do not
    # compare equal. exact=True specifies that the type of the indices must be
    # exactly the same
    pdt.assert_index_equal(left.index, right.index,
                           exact=True,
                           check_names=True)
