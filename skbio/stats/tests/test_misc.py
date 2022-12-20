# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

from skbio.stats._misc import _pprint_strs


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
