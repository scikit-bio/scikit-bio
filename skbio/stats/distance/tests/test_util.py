# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from skbio.stats.distance._util import is_id_pair


class TestIsIDPair(unittest.TestCase):
    def test_positive(self):
        self.assertTrue(is_id_pair(('', '')))
        self.assertTrue(is_id_pair(('', '123')))
        self.assertTrue(is_id_pair(('a', 'b')))
        self.assertTrue(is_id_pair(('xyz', 'ab')))
        self.assertTrue(is_id_pair(('xyz', 'xyz')))

    def test_negative_wrong_container_type(self):
        self.assertFalse(is_id_pair(['abc', 'xyz']))
        self.assertFalse(is_id_pair({'abc', 'xyz'}))
        self.assertFalse(is_id_pair('abc'))

    def test_negative_wrong_element_type(self):
        self.assertFalse(is_id_pair(('xyz', 123)))
        self.assertFalse(is_id_pair((123, 'xyz')))
        self.assertFalse(is_id_pair((123, 123)))
        self.assertFalse(is_id_pair(('abc', 42.5)))
        self.assertFalse(is_id_pair(('abc', ['xyz'])))

    def test_negative_wrong_length(self):
        self.assertFalse(is_id_pair(()))
        self.assertFalse(is_id_pair(('a',)))
        self.assertFalse(is_id_pair(('a', 'b', 'c')))


if __name__ == '__main__':
    unittest.main()
