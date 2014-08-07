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


if __name__ == '__main__':
    main()
