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


class IDValidationTests(object):
    """Validate IDs in scikit-bio objects.

    Must define a 'id_cls' instance variable in setUp.
    Must define extra 'id_kwargs', and potentially an 'id_args' if the object takes
    an undefined number of positional arguments.

    """

    def create_instance(self, id=None):
        if id is None:
            instance = self.id_cls(**self.id_kwargs)
        else:
            instance = self.id_cls(id=id, **self.id_kwargs)

        return instance

    def test_no_id_set(self):
        expected_id = ''
        instance = self.create_instance()
        self.assertEqual(instance.id, expected_id)

    def test_valid_id_construction_and_access(self):
        expected_id = 'whatever'
        instance = self.create_instance(expected_id)
        self.assertEqual(instance.id, expected_id)

    def test_invalid_id_non_string_construction(self):
        invalid_id = ('a', 1)

        with self.assertRaises(ValueError) as e:
            self.create_instance(invalid_id)
        self.assertIn('ID', str(e.exception))

    def test_set_id(self):
        valid_id = 'id'

        with self.assertRaises(AttributeError):
            instance = self.create_instance()
            instance.id = valid_id
