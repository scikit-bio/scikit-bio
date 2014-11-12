# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from unittest import TestCase, main

import numpy as np

from skbio.stats import p_value_to_str
from skbio.stats._misc import _pprint_strs


class PValueToStrTests(TestCase):
    def setUp(self):
        self.p_value = 0.119123123123

    def test_valid_input(self):
        obs = p_value_to_str(self.p_value, 100)
        self.assertEqual(obs, '0.12')

        obs = p_value_to_str(self.p_value, 250)
        self.assertEqual(obs, '0.12')

        obs = p_value_to_str(self.p_value, 1000)
        self.assertEqual(obs, '0.119')

        obs = p_value_to_str(0.0055623489, 999)
        self.assertEqual(obs, '0.006')

    def test_too_few_permutations(self):
        obs = p_value_to_str(self.p_value, 9)
        self.assertEqual(obs, 'Too few permutations to compute p-value '
                              '(permutations = 9)')

        obs = p_value_to_str(self.p_value, 1)
        self.assertEqual(obs, 'Too few permutations to compute p-value '
                              '(permutations = 1)')

        obs = p_value_to_str(self.p_value, 0)
        self.assertEqual(obs, 'Too few permutations to compute p-value '
                              '(permutations = 0)')

    def test_missing_or_invalid_p_value(self):
        obs = p_value_to_str(None, 0)
        self.assertEqual(obs, 'N/A')

        obs = p_value_to_str(np.nan, 0)
        self.assertEqual(obs, 'N/A')


class PPrintStrsTests(TestCase):
    def test_truncation(self):
        # truncation between items (on comma)
        exp = "'a', ..."
        obs = _pprint_strs(['a', 'b', 'c'], max_chars=4)
        self.assertEqual(obs, exp)

        # truncation between items (on space)
        exp = "'a', ..."
        obs = _pprint_strs(['a', 'b', 'c'], max_chars=5)
        self.assertEqual(obs, exp)

        # truncation on item
        exp = "'a', ..."
        obs = _pprint_strs(['a', 'b', 'c'], max_chars=6)
        self.assertEqual(obs, exp)

        # truncation (no items)
        exp = "..."
        obs = _pprint_strs(['a', 'b', 'c'], max_chars=2)
        self.assertEqual(obs, exp)

    def test_no_truncation(self):
        exp = "'a'"
        obs = _pprint_strs(['a'], max_chars=3)
        self.assertEqual(obs, exp)

        exp = "'a', 'b', 'c'"
        obs = _pprint_strs(['a', 'b', 'c'])
        self.assertEqual(obs, exp)

        exp = "'a', 'b', 'c'"
        obs = _pprint_strs(['a', 'b', 'c'], max_chars=13)
        self.assertEqual(obs, exp)

    def test_non_default_delimiter_and_suffix(self):
        exp = "'abc','defg',...."
        obs = _pprint_strs(['abc', 'defg', 'hi', 'jklmno'], max_chars=14,
                           delimiter=',', suffix='....')
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
