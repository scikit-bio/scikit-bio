from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from os import remove

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
        pass

    def tearDown(self):
        pass

    def retrieval_test(self):
        class Google(URLGetter):
            base_url = 'http://www.google.com'
        g = Google()

        # test URL construction
        self.assertEqual(str(g), 'http://www.google.com')
        # test reading
        text = g.read()
        self.assertTrue('<title>Google</title>' in text)
        # test file getting
        fname = '/tmp/google_test'
        g.retrieve(fname)
        g_file = open(fname)
        g_text = g_file.read()
        self.assertEqual(g_text, text)
        g_text.close()
        remove(fname)

if __name__ == '__main__':
    main()
