#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main

from skbio.util.sort import signed_natsort, natsort


class SortTests(TestCase):
    """Test object for sort functions"""

    def setUp(self):
        """"""
        pass

    def test_natsort(self):
        """natsort should perform numeric comparisons on strings"""
        # string with alpha and numerics sort correctly
        s = ['sample1', 'sample2', 'sample11', 'sample12']
        exp = ['sample1', 'sample2', 'sample11', 'sample12']
        self.assertEqual(natsort(s), exp)
        s.reverse()
        self.assertEqual(natsort(s), exp)
        self.assertEqual(natsort(list('cba321')), list('123abc'))

        # strings with alpha only sort correctly
        self.assertEqual(natsort(list('cdba')), list('abcd'))

        # string of ints sort correctly
        s = ['11', '2', '1', '0']
        exp = ['0', '1', '2', '11']
        self.assertEqual(natsort(s), exp)

        # strings of floats sort correctly
        s = ['1.11', '1.12', '1.00', '0.009']
        exp = ['0.009', '1.00', '1.11', '1.12']
        self.assertEqual(natsort(s), exp)

        # tuples sort correctly
        s = [('11', 'A'), ('2', 'B'), ('1', 'C'), ('0', 'D')]
        exp = [('0', 'D'), ('1', 'C'), ('2', 'B'), ('11', 'A')]
        self.assertEqual(natsort(s), exp)

    def test_natsort_case_insensitive(self):
        """natsort should perform numeric comparisons on strings and is _not_
        case-sensitive"""

        # string with alpha and numerics sort correctly
        s = ['sample1', 'sample2', 'sample11', 'sample12',
             'SAmple1', 'Sample2']

        # expected values
        exp_natsort = ['SAmple1', 'Sample2', 'sample1', 'sample2', 'sample11',
                       'sample12']
        exp_natsort_case_insensitive = ['sample1', 'SAmple1', 'sample2',
                                        'Sample2', 'sample11', 'sample12']

        # test natsort
        self.assertEqual(natsort(s), exp_natsort)
        # test natsort case insensitive
        self.assertEqual(natsort(s, case_sensitive=False),
                         exp_natsort_case_insensitive)
        s.reverse()
        # test natsort
        self.assertEqual(natsort(s), exp_natsort)
        # test natsort case insensitive
        exp_natsort_case_insensitive = ['SAmple1', 'sample1', 'Sample2',
                                        'sample2', 'sample11', 'sample12']
        self.assertEqual(natsort(s, case_sensitive=False),
                         exp_natsort_case_insensitive)
        self.assertEqual(natsort(list('cbaA321')), list('123Aabc'))
        self.assertEqual(natsort(list('cbaA321'), case_sensitive=False),
                         list('123aAbc'))

        # strings with alpha only sort correctly
        self.assertEqual(natsort(list('cdBa'), case_sensitive=False),
                         list('aBcd'))

        # string of ints sort correctly
        s = ['11', '2', '1', '0']
        exp = ['0', '1', '2', '11']
        self.assertEqual(natsort(s, case_sensitive=False), exp)

        # strings of floats sort correctly
        s = ['1.11', '1.12', '1.00', '0.009']
        exp = ['0.009', '1.00', '1.11', '1.12']
        self.assertEqual(natsort(s, case_sensitive=False), exp)

        # tuples sort correctly
        s = [('11', 'A'), ('2', 'B'), ('1', 'C'), ('0', 'D')]
        exp = [('0', 'D'), ('1', 'C'), ('2', 'B'), ('11', 'A')]
        self.assertEqual(natsort(s, case_sensitive=False), exp)

    def test_signed_sort(self):
        """Test correct sorting of different data types"""

        # an empty list must be returned when an empty list needs to be sorted
        self.assertEqual(signed_natsort([]), [])

        # tuples that can be sorted by type-casting the first element
        test_list = [('9', 'SampleA'), ('-1', 'SampleD'), ('7', 'SampleC'),
                     ('-2', 'SampleE'), ('-0.11', 'SampleF'),
                     ('17.11', 'SampleB'), ('100', 'SampleG'),
                     ('13', 'SampleH')]
        expected_result = [('-2', 'SampleE'), ('-1', 'SampleD'),
                           ('-0.11', 'SampleF'), ('7', 'SampleC'),
                           ('9', 'SampleA'), ('13', 'SampleH'),
                           ('17.11', 'SampleB'), ('100', 'SampleG')]
        output = signed_natsort(test_list)
        self.assertEquals(output, expected_result)

        # tuples that must be sorted alphabetically
        test_list = [('Cygnus', 'SampleA'), ('Cepheus', 'SampleD'),
                     ('Auriga', 'SampleC'), ('Grus', 'SampleE'),
                     ('Hydra', 'SampleF'), ('Carina', 'SampleB'),
                     ('Orion', 'SampleG'), ('Lynx', 'SampleH')]
        expected_result = [('Auriga', 'SampleC'), ('Carina', 'SampleB'),
                           ('Cepheus', 'SampleD'), ('Cygnus', 'SampleA'),
                           ('Grus', 'SampleE'), ('Hydra', 'SampleF'),
                           ('Lynx', 'SampleH'), ('Orion', 'SampleG')]
        output = signed_natsort(test_list)
        self.assertEquals(output, expected_result)

        # mixed case, tuples will be sorted alpha-numerically
        test_list = [('Cygnus', 'SampleA'), ('Cepheus', 'SampleD'),
                     ('Auriga', 'SampleC'), ('Grus', 'SampleE'),
                     ('-0.11', 'SampleF'), ('17.11', 'SampleB'),
                     ('100', 'SampleG'), ('Lynx', 'SampleH')]
        expected_result = [('17.11', 'SampleB'), ('100', 'SampleG'),
                           ('-0.11', 'SampleF'), ('Auriga', 'SampleC'),
                           ('Cepheus', 'SampleD'), ('Cygnus', 'SampleA'),
                           ('Grus', 'SampleE'), ('Lynx', 'SampleH')]
        output = signed_natsort(test_list)
        self.assertEquals(output, expected_result)

        # mixed case just a list
        test_list = ['foo', 'bar', '-100', '12', 'spam', '4', '-1']
        expected_result = ['4', '12', '-1', '-100', 'bar', 'foo', 'spam']
        output = signed_natsort(test_list)
        self.assertEquals(output, expected_result)

        # list of elements that can be type-casted
        test_list = ['0', '1', '14', '12', '-15', '4', '-1']
        expected_result = ['-15', '-1', '0', '1', '4', '12', '14']
        output = signed_natsort(test_list)
        self.assertEquals(output, expected_result)

        # mixed dict case
        test_dict = {'foo': 'a', 'bar': 'b', '-100': '1', '12': '11',
                     'spam': 'q', '4': '11', '-1': 'e'}
        expected_result = ['4', '12', '-1', '-100', 'bar', 'foo', 'spam']
        output = signed_natsort(test_dict)
        self.assertEquals(output, expected_result)

        # dict where the keys can be type-casted
        test_dict = {'0': 'foo', '1': 'bar', '14': 'stand', '12': 'eggs',
                     '-15': 'q', '4': 'b', '-1': 'h'}
        expected_result = ['-15', '-1', '0', '1', '4', '12', '14']
        output = signed_natsort(test_dict)
        self.assertEquals(output, expected_result)


if __name__ == '__main__':
    main()
