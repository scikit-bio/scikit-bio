# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import unittest

from skbio.io._base import _chunk_str, _int_to_ordinal_str


class ChunkStrTests(unittest.TestCase):
    def test_even_split(self):
        self.assertEqual(_chunk_str('abcdef', 6, ' '), 'abcdef')
        self.assertEqual(_chunk_str('abcdef', 3, ' '), 'abc def')
        self.assertEqual(_chunk_str('abcdef', 2, ' '), 'ab cd ef')
        self.assertEqual(_chunk_str('abcdef', 1, ' '), 'a b c d e f')
        self.assertEqual(_chunk_str('a', 1, ' '), 'a')
        self.assertEqual(_chunk_str('abcdef', 2, ''), 'abcdef')

    def test_no_split(self):
        self.assertEqual(_chunk_str('', 2, '\n'), '')
        self.assertEqual(_chunk_str('a', 100, '\n'), 'a')
        self.assertEqual(_chunk_str('abcdef', 42, '|'), 'abcdef')

    def test_uneven_split(self):
        self.assertEqual(_chunk_str('abcdef', 5, '|'), 'abcde|f')
        self.assertEqual(_chunk_str('abcdef', 4, '|'), 'abcd|ef')
        self.assertEqual(_chunk_str('abcdefg', 3, ' - '), 'abc - def - g')

    def test_invalid_n(self):
        with self.assertRaisesRegexp(ValueError, 'n=0'):
            _chunk_str('abcdef', 0, ' ')

        with self.assertRaisesRegexp(ValueError, 'n=-42'):
            _chunk_str('abcdef', -42, ' ')


class IntToOrdinalStrTests(unittest.TestCase):
    def test_valid_range(self):
        # taken and modified from http://stackoverflow.com/a/20007730/3776794
        exp = ['0th', '1st', '2nd', '3rd', '4th', '5th', '6th', '7th', '8th',
               '9th', '10th', '11th', '12th', '13th', '14th', '15th', '16th',
               '17th', '18th', '19th', '20th', '21st', '22nd', '23rd', '24th',
               '25th', '26th', '27th', '28th', '29th', '30th', '31st', '32nd',
               '100th', '101st']
        obs = [_int_to_ordinal_str(n) for n in range(0, 33) + [100, 101]]
        self.assertEqual(obs, exp)

    def test_invalid_n(self):
        with self.assertRaisesRegexp(ValueError, '-1'):
            _int_to_ordinal_str(-1)


if __name__ == '__main__':
    unittest.main()
