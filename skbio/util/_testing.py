# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import inspect

import pandas.util.testing as pdt
from nose import core
from nose.tools import nottest
from future.utils import PY3

from ._decorator import experimental


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
        if not PY3:
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
