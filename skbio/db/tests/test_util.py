#!/usr/bin/env python
"""Tests of the db utility functions and classes."""
from unittest import TestCase, main
from skbio.db.util import UrlGetter, expand_slice, last_nondigit_index,make_lists_of_expanded_slices_of_set_size,make_lists_of_accessions_of_set_size
from os import remove

class db_util_tests(TestCase):
    """Tests of top-level functions."""
    def test_last_nondigit_index(self):
        """last_nondigit_index should return i such that s[i:] is numeric"""
        ldi = last_nondigit_index
        self.assertEqual(ldi('3'), 0)
        self.assertEqual(ldi('345'), 0)
        self.assertEqual(ldi('a34'), 1)
        self.assertEqual(ldi('abc234'), 3)
        self.assertEqual(ldi('abcd'), None)

    def test_expand_slice(self):
        """expand_slice should get accession range"""
        self.assertEqual(expand_slice(slice('AF1001','AF1003')), \
            ['AF1001','AF1002','AF1003'])
        #can't expand if accession prefixes
        self.assertRaises(TypeError, expand_slice, slice('AF100:','AG1002'))
        #should keep leading zeros
        self.assertEqual(expand_slice(slice('AF0001','AF0003')), \
            ['AF0001','AF0002','AF0003'])

    def test_make_lists_of_expanded_slices_of_set_size(self):
        """make_lists_of_expanded_slices_of_set_size: should return a
        list of lists"""
        expected_list = ['HM780503 HM780504 HM780505','HM780506']
        observed = make_lists_of_expanded_slices_of_set_size(slice('HM780503','HM780506'),size_limit=3)
        self.assertEqual(observed,expected_list)

    def make_lists_of_accessions_of_set_size(self):
        """make_lists_of_expanded_slices_of_set_size: should return a
        list of lists"""
        expected_list = ['HM780503 HM780506 HM780660 HM780780']
        observed = make_lists_of_accessions_of_set_size(['HM780503','HM780506', 'HM780660', 'HM780780'],size_limit=3)        
        self.assertEqual(observed,expected_list)


class UrlGetterTests(TestCase):
    """Tests of the UrlGetter class"""
    def retrieval_test(self):
        """Urlgetter should init, read and retrieve"""
        class Google(UrlGetter):
            BaseUrl='http://www.google.com'
        g = Google()
        #test URL construction
        self.assertEqual(str(g), 'http://www.google.com')
        #test reading
        text = g.read()
        assert '<title>Google</title>' in text
        #test file getting
        fname = '/tmp/google_test'
        g.retrieve(fname)
        g_file = open(fname)
        g_text = g_file.read()
        self.assertEqual(g_text, text)
        g_text.close()
        remove(fname)

if __name__ == '__main__':
    main()
