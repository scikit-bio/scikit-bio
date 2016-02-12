# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import unittest

from skbio.sequence import IntervalMetadataMixin, IntervalMetadata, Interval, Feature
from skbio.sequence._interval import _polish_interval


class TestIntervalMetadataMixin(unittest.TestCase):
    def setUp(self):
        pass

    def test_init(self):
        d = {Feature(gene='sagA'): [(0, 2), (4, 7)],
             Feature(gene='sagB'): [(3, 5)]}
        md = IntervalMetadata(features=d)
        self.assertEquals(md.features, d)

    def test_polish_interval_empty(self):
        res = _polish_interval(())
        self.assertTrue(res is None)

    def test_polish_interval_tuple(self):
        st, end = _polish_interval((1, 2))
        self.assertEquals(st, 1)
        self.assertEquals(end, 2)

    def test_polish_interval_interval(self):
        st, end = _polish_interval(Interval(1, 2))
        self.assertEquals(st, 1)
        self.assertEquals(end, 2)

    def test_polish_interval_point(self):
        st, end = _polish_interval(1)
        self.assertEquals(st, 1)
        self.assertEquals(end, 2)

    def test_query(self):
        interval_metadata = IntervalMetadata(features={
                           Feature(gene='sagA'): [(0, 2), (4, 7)],
                           Feature(gene='sagB'): [(3, 5)]
                       })

        feats = interval_metadata.query((1, 2))
        self.assertEqual(feats, [Feature(gene='sagA')])

        feats = interval_metadata.query(gene='sagB')
        self.assertEqual(feats, [Feature(gene='sagB')])

    def test_query_duplicate1(self):
        interval_metadata = IntervalMetadata(features={
                           Feature(gene='sagA'): [(0, 2), (4, 7)],
                           Feature(gene='sagB'): [(3, 5)]
                       })

        feats = interval_metadata.query((1, 5))
        self.assertEqual(set(feats),
                         {Feature(gene='sagA'),
                          Feature(gene='sagB')})

    def test_query_duplicate2(self):
        interval_metadata = IntervalMetadata(features={
                           Feature(gene='sagA'): [(0, 2), (4, 7)],
                           Feature(gene='sagB'): [(3, 5)],
                           Feature(gene='sagB'): [(14, 29)]
                       })
        feats = interval_metadata.query(gene='sagB')
        self.assertEqual(feats, [Feature(gene='sagB')])


    def test_add(self):
        interval_metadata = IntervalMetadata()
        interval_metadata.add(Feature(gene='sagA'), 1, (4, 7))
        interval_metadata.add(Feature(gene='sagB'), (3, 5))

        # Relies on the test_query method to work
        feats = interval_metadata.query((1, 2))
        self.assertEqual(feats, [Feature(gene='sagA')])

        feats = interval_metadata.query(gene='sagB')
        self.assertEqual(feats, [Feature(gene='sagB')])

    def test_add_empty(self):
        interval_metadata = IntervalMetadata()
        interval_metadata.add(Feature(gene='sagA'))
        feats = interval_metadata.query(gene='sagA')
        self.assertEqual(feats, [Feature(gene='sagA')])

        interval_metadata.add(Feature(gene='sagB'), ())
        feats = interval_metadata.query(gene='sagB')
        self.assertEqual(feats, [Feature(gene='sagB')])

if __name__=='__main__':
    unittest.main()
