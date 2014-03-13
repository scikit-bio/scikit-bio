#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, bipy development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from StringIO import StringIO
from tempfile import NamedTemporaryFile
from os.path import exists

from bipy.util.misc import safe_md5, remove_files
from bipy.util.unit_test import TestCase, main


class UtilTests(TestCase):
    """Test object for the miscellaneous utility functions"""
    def test_safe_md5(self):
        """Make sure we have the expected md5"""
        exp = 'ab07acbb1e496801937adfa772424bf7'

        fd = StringIO('foo bar baz')
        obs = safe_md5(fd)
        self.assertEqual(obs.hexdigest(), exp)

        fd.close()

    def test_remove_files(self):
        """Remove files functions as expected """
        # create list of temp file paths
        test_fds = [NamedTemporaryFile(delete=False) for i in range(5)]
        test_filepaths = [element.name for element in test_fds]

        # should work just fine
        remove_files(test_filepaths)

        # check that an error is raised on trying to remove the files...
        self.assertRaises(OSError, remove_files, test_filepaths)

        # touch one of the filepaths so it exists
        extra_file = NamedTemporaryFile(delete=False).name
        test_filepaths.append(extra_file)

        # no error is raised on trying to remove the files 
        # (although 5 don't exist)...
        remove_files(test_filepaths, error_on_missing=False)
        # ... and the existing file was removed
        self.assertFalse(exists(extra_file))

        # try to remove them with remove_files and verify that an IOError is 
        # raises
        self.assertRaises(OSError, remove_files, test_filepaths)

        # now get no error when error_on_missing=False
        remove_files(test_filepaths,error_on_missing=False)

if __name__ == '__main__':
    main()
