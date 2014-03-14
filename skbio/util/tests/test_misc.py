#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from StringIO import StringIO
from tempfile import NamedTemporaryFile
from os.path import exists, join
from skbio.util.misc import safe_md5, remove_files
from unittest import TestCase, main
from shutil import rmtree
from tempfile import mkdtemp
from uuid import uuid4

from skbio.util.misc import (safe_md5, remove_files, create_dir)


class MiscTests(TestCase):

    """Test object for the miscellaneous utility functions"""

    def setUp(self):
        self.dirs_to_remove = []

    def tearDown(self):
        for element in self.dirs_to_remove:
            rmtree(element)

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
        remove_files(test_filepaths, error_on_missing=False)

    def test_create_dir(self):
        """create_dir creates dir and fails meaningful."""

        # create a directory
        tmp_dir_path = mkdtemp()

        # create a random temporary directory name
        tmp_dir_path2 = join(mkdtemp(), str(uuid4()))
        tmp_dir_path3 = join(mkdtemp(), str(uuid4()))

        self.dirs_to_remove += [tmp_dir_path, tmp_dir_path2, tmp_dir_path3]

        # create on existing dir raises OSError if fail_on_exist=True
        self.assertRaises(OSError, create_dir, tmp_dir_path,
                          fail_on_exist=True)
        self.assertEquals(create_dir(tmp_dir_path, fail_on_exist=True,
                                     handle_errors_externally=True), 1)

        # return should be 1 if dir exist and fail_on_exist=False
        self.assertEqual(create_dir(tmp_dir_path, fail_on_exist=False), 1)

        # if dir not there make it and return always 0
        self.assertEqual(create_dir(tmp_dir_path2), 0)
        self.assertEqual(create_dir(tmp_dir_path3, fail_on_exist=True), 0)

if __name__ == '__main__':
    main()
