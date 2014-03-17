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
from os import rmdir
from os.path import exists, join
from unittest import TestCase, main
from shutil import rmtree
from tempfile import mkdtemp
from uuid import uuid4

from skbio.util.misc import (safe_md5, remove_files, create_dir,
                             get_random_directory_name)


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

    def test_get_random_directory_name(self):
        """get_random_directory_name functions as expected """
        # repeated calls yield different directory names
        dirs = []
        for i in range(100):
            d = get_random_directory_name(suppress_mkdir=True)
            self.assertTrue(d not in dirs)
            dirs.append(d)

        actual = get_random_directory_name(suppress_mkdir=True)
        self.assertFalse(exists(actual),'Random dir exists: %s' % actual)
        self.assertTrue(actual.startswith('/'),\
         'Random dir is not a full path: %s' % actual)

        # prefix, suffix and output_dir are used as expected
        actual = get_random_directory_name(suppress_mkdir=True,prefix='blah',\
            output_dir='/tmp/',suffix='stuff')
        self.assertTrue(actual.startswith('/tmp/blah2'),\
         'Random dir does not begin with output_dir + prefix '+\
         '+ 2 (where 2 indicates the millenium in the timestamp): %s' % actual)
        self.assertTrue(actual.endswith('stuff'),\
         'Random dir does not end with suffix: %s' % actual)

        # changing rand_length functions as expected
        actual1 = get_random_directory_name(suppress_mkdir=True)
        actual2 = get_random_directory_name(suppress_mkdir=True,\
            rand_length=10)
        actual3 = get_random_directory_name(suppress_mkdir=True,\
            rand_length=0)
        self.assertTrue(len(actual1) > len(actual2) > len(actual3),\
         "rand_length does not affect directory name lengths "+\
         "as expected:\n%s\n%s\n%s" % (actual1,actual2,actual3))

        # changing the timestamp pattern functions as expected
        actual1 = get_random_directory_name(suppress_mkdir=True)
        actual2 = get_random_directory_name(suppress_mkdir=True,\
            timestamp_pattern='%Y')
        self.assertNotEqual(actual1,actual2)
        self.assertTrue(len(actual1)>len(actual2),\
            'Changing timestamp_pattern does not affect directory name')
        # empty string as timestamp works
        actual3 = get_random_directory_name(suppress_mkdir=True,\
            timestamp_pattern='')
        self.assertTrue(len(actual2) > len(actual3))

        # creating the directory works as expected
        actual = get_random_directory_name(output_dir='/tmp/',\
         prefix='get_random_directory_test')
        self.assertTrue(exists(actual))
        rmdir(actual)

if __name__ == '__main__':
    main()
