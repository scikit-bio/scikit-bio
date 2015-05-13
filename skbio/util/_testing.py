# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import inspect
from nose import core
from nose.tools import nottest
from future.utils import PY3
import numpy.testing as npt
from pandas.util.testing import assert_frame_equal, assert_series_equal


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


def assert_ordination_results_equal(left, right, ignore_method_name=False):
    """Assert that ordination results objects are equal.

    This is a helper function intended to be used in unit tests that need to
    compare ``OrdinationResults`` objects.

    For numeric attributes (e.g., eigvals, site, etc.),
    ``numpy.testing.assert_almost_equal`` is used. Otherwise,
    ``numpy.testing.assert_equal`` is used for comparisons. An assertion is
    in place to ensure the two objects are exactly the same type.

    Parameters
    ----------
    left, right : OrdinationResults
        Ordination results to be compared for equality.

    Raises
    ------
    AssertionError
        If the two objects are not equal.

    """
    npt.assert_equal(type(left) is type(right), True)

    if not ignore_method_name:
        npt.assert_equal(left.short_method_name, right.short_method_name)
        npt.assert_equal(left.long_method_name, right.long_method_name)

    # eigvals should always be present
    assert_series_equal(left.eigvals, right.eigvals)

    # samples should always be present
    assert_frame_equal(left.samples, right.samples)

    # assert_frame_equal doesn't like None...
    def _assert_optional_frame_equal(left, right):
        if left is None or right is None:
            npt.assert_equal(left, right)
        else:
            assert_frame_equal(left, right)

    _assert_optional_frame_equal(left.features, right.features)
    _assert_optional_frame_equal(left.biplot_scores, right.biplot_scores)
    _assert_optional_frame_equal(left.sample_constraints,
                                 right.sample_constraints)

    # assert_series_equal doesn't like None...
    def _assert_optional_series_equal(left, right):
        if left is None or right is None:
            npt.assert_equal(left, right)
        else:
            assert_series_equal(left, right)

    _assert_optional_series_equal(left.proportion_explained,
                                  right.proportion_explained)
