# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import unittest

from skbio.metadata import Feature


class TestFeature(unittest.TestCase):
    def setUp(self):
        self.exp = [
            {},
            {'db_xref': '"taxon:562"',
             'index_': 0,
             'left_partial_': False},
            {'codon_start': '1',
             'db_xref': ('"GI:145230"', '"taxon:562"'),
             'index_': 1}]

    def test_repr(self):
        for i in self.exp:
            exp = ';'.join('{0}:{1}'.format(k, i[k]) for k in i)
            self.assertEqual(repr(Feature(**i)), exp)

    def test_getitem(self):
        for i in self.exp:
            x = Feature(**i)
            for j in x:
                self.assertEqual(i[j], x[j])

    def test_hash(self):
        for i in self.exp:
            x1 = Feature(**i)
            x2 = Feature(**i)
            self.assertEqual(hash(x1), hash(x2))
            i['foo'] = 'spam'
            x3 = Feature(**i)
            self.assertNotEqual(hash(x1), hash(x3))

    def test_update(self):
        for i in self.exp:
            x1 = Feature(**i)
            x2 = Feature(**i)

            i['k'] = 10
            exp = Feature(**i)
            x1 = x1.update(k=10)
            self.assertEqual(hash(x1), hash(exp))


    ## TODO: Need tests for keys(), items(), values()

if __name__ == '__main__':
    unittest.main()
