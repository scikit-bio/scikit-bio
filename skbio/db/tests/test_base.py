from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from tempfile import mkstemp
from os import remove, close

from skbio.db.base import (URLGetter, expand_slice, last_nondigit_index,
                           make_lists_of_expanded_slices_of_set_size,
                           make_lists_of_accessions_of_set_size)


class BaseTests(TestCase):
    def test_last_nondigit_index(self):
        ldi = last_nondigit_index
        self.assertEqual(ldi('3'), 0)
        self.assertEqual(ldi('345'), 0)
        self.assertEqual(ldi('a34'), 1)
        self.assertEqual(ldi('abc234'), 3)
        self.assertEqual(ldi('abcd'), None)

    def test_expand_slice(self):
        self.assertEqual(expand_slice(slice('AF1001', 'AF1003')),
                         ['AF1001', 'AF1002', 'AF1003'])

        # can't expand if accession prefixes
        self.assertRaises(TypeError, expand_slice, slice('AF100:', 'AG1002'))

        # should keep leading zeros
        self.assertEqual(expand_slice(slice('AF0001', 'AF0003')),
                         ['AF0001', 'AF0002', 'AF0003'])

    def test_make_lists_of_expanded_slices_of_set_size(self):
        expected_list = ['HM780503 HM780504 HM780505', 'HM780506']
        observed = make_lists_of_expanded_slices_of_set_size(
            slice('HM780503', 'HM780506'), size_limit=3)
        self.assertEqual(observed, expected_list)

    def make_lists_of_accessions_of_set_size(self):
        expected_list = ['HM780503 HM780506 HM780660 HM780780']
        observed = make_lists_of_accessions_of_set_size(
            ['HM780503', 'HM780506', 'HM780660', 'HM780780'], size_limit=3)
        self.assertEqual(observed, expected_list)


class URLGetterTests(TestCase):
    def setUp(self):
        fd, self.input_fp = mkstemp(prefix="google_file_", suffix='.temp')
        close(fd)

        self.g = URLGetter()
        self.g.base_url = 'http://www.google.com'

    def tearDown(self):
        try:
            remove(self.input_fp)
        except OSError:
            pass

    def test_str(self):
        # test URL construction
        self.assertEqual(str(self.g), 'http://www.google.com')

    def test_read(self):
        # test reading
        text = self.g.read()
        self.assertTrue('<title>Google</title>' in str(text))

    def test_retrieve(self):
        # test file getting
        self.g.retrieve(self.input_fp)
        g_file = open(self.input_fp)
        g_text = g_file.read()
        self.assertTrue('<title>Google</title>' in g_text)
        g_file.close()

if __name__ == '__main__':
    main()
