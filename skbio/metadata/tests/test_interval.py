# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import unittest

from skbio.metadata import IntervalMetadata, Interval, Feature
from skbio.metadata._interval import _polish_interval


class TestIntervalMetadataMixin(unittest.TestCase):
    def setUp(self):
        pass

    def test_init(self):
        d = {Feature(gene='sagA', location=0): [(0, 2), (4, 7)],
             Feature(gene='sagB', location=0): [(3, 5)]}
        md = IntervalMetadata(features=d)
        self.assertEquals(md.features, d)

    def test_polish_interval_empty(self):
        res = _polish_interval(())
        self.assertTrue(res is None)

    def test_polish_interval_tuple(self):
        st, end = _polish_interval((1, 2))
        self.assertEqual(st, 1)
        self.assertEqual(end, 2)

    def test_polish_interval_interval(self):
        st, end = _polish_interval(Interval(1, 2))
        self.assertEqual(st, 1)
        self.assertEqual(end, 2)

    def test_polish_interval_point(self):
        st, end = _polish_interval(1)
        self.assertEqual(st, 1)
        self.assertEqual(end, 2)

    def test_query(self):
        interval_metadata = IntervalMetadata(features={
                           Feature(gene='sagA', location=0): [(0, 2), (4, 7)],
                           Feature(gene='sagB', location=0): [(3, 5)]
                       })

        feats = interval_metadata.query((1, 2))
        self.assertEqual(feats, [Feature(gene='sagA', location=0)])

        feats = interval_metadata.query(gene='sagB')
        self.assertEqual(feats, [Feature(gene='sagB', location=0)])

    def test_query_duplicate1(self):
        interval_metadata = IntervalMetadata(features={
                           Feature(gene='sagA', location=0): [(0, 2), (4, 7)],
                           Feature(gene='sagB', location=0): [(3, 5)]
                       })

        feats = interval_metadata.query((1, 5))
        self.assertEqual(set(feats),
                         {Feature(gene='sagA', location=0),
                          Feature(gene='sagB', location=0)})

    def test_query_duplicate2(self):
        interval_metadata = IntervalMetadata(features={
                           Feature(gene='sagA', location=0): [(0, 2), (4, 7)],
                           Feature(gene='sagB', location=0): [(3, 5)],
                           Feature(gene='sagB', location=0): [(14, 29)]
                       })
        feats = interval_metadata.query(gene='sagB')
        self.assertEqual(feats, [Feature(gene='sagB', location=0)])

    def test_reverse_complement(self):
        interval_metadata = IntervalMetadata()
        interval_metadata.add(Feature(gene='sagB', location=0), (3, 5))
        iv = interval_metadata.reverse_complement(length=10)
        feats = iv.query((5, 7))
        self.assertEqual(feats, [Feature(gene='sagB', location=0)])

    def test_update(self):
        interval_metadata = IntervalMetadata()
        interval_metadata.add(Feature(gene='sagA', location=0), 1, (4, 7))
        interval_metadata.add(Feature(gene='sagB', location=0), (3, 5))
        interval_metadata.update(Feature(gene='sagB', location=0),
                                 Feature(gene='sagB', location=1))
        feats = interval_metadata.query(gene='sagB')
        self.assertEqual(feats, [Feature(gene='sagB', location=1)])

    def test_update_multiple(self):
        # Tests to see if all corresponding intervals are updated
        interval_metadata = IntervalMetadata()
        interval_metadata.add(Feature(gene='sagA', location=0),
                              1, (4, 7))
        interval_metadata.add(Feature(gene='sagB', location=0),
                              (3, 5), (10, 15))
        interval_metadata.update(Feature(gene='sagB', location=0),
                                 Feature(gene='sagB', location=1))
        feats = interval_metadata.query(gene='sagB')
        self.assertEqual(feats, [Feature(gene='sagB', location=1)])

        feats = sorted(interval_metadata.query((3, 5)),
                       key=lambda x: x['location'])
        self.assertEqual(feats, sorted([Feature(gene='sagA', location=0),
                                        Feature(gene='sagB', location=1)],
                                       key=lambda x: x['location']))

        feats = interval_metadata.query((10, 15))
        self.assertEqual(feats, [Feature(gene='sagB', location=1)])

    def test_add(self):
        interval_metadata = IntervalMetadata()
        interval_metadata.add(Feature(gene='sagA', location=0), 1, (4, 7))
        interval_metadata.add(Feature(gene='sagB', location=0), (3, 5))

        # Relies on the test_query method to work
        feats = interval_metadata.query((1, 2))
        self.assertEqual(feats, [Feature(gene='sagA', location=0)])

        feats = interval_metadata.query(gene='sagB')
        self.assertEqual(feats, [Feature(gene='sagB', location=0)])

    def test_add_empty(self):
        interval_metadata = IntervalMetadata()
        interval_metadata.add(Feature(gene='sagA', location=0))
        feats = interval_metadata.query(gene='sagA')
        self.assertEqual(feats, [Feature(gene='sagA', location=0)])

        interval_metadata.add(Feature(gene='sagB', location=0), ())
        feats = interval_metadata.query(gene='sagB')
        self.assertEqual(feats, [Feature(gene='sagB', location=0)])

    def test_concat(self):
        interval_metadata1 = IntervalMetadata(features={
            Feature(gene='sagA', location='0'): [(0, 2), (4, 7)],
            Feature(gene='sagB', location='3'): [(3, 5)]
        })

        interval_metadata2 = IntervalMetadata(features={
                   Feature(gene='sagC', location='10'): [(10, 12), (24, 27)],
                   Feature(gene='sagD', location='13'): [(13, 15)]
               })
        catted_metadata = interval_metadata1.concat(interval_metadata2)
        self.assertEqual(catted_metadata.features, {
                         Feature(gene='sagA', location='0'): [(0, 2), (4, 7)],
                         Feature(gene='sagB', location='3'): [(3, 5)],
                         Feature(gene='sagC', location='10'): [(10, 12),
                                                               (24, 27)],
                         Feature(gene='sagD', location='13'): [(13, 15)]})

    def test_concat_inplace(self):
        interval_metadata1 = IntervalMetadata(features={
                   Feature(gene='sagA', location='0'): [(0, 2), (4, 7)],
                   Feature(gene='sagB', location='3'): [(3, 5)]
               })

        interval_metadata2 = IntervalMetadata(features={
                   Feature(gene='sagC', location='10'): [(10, 12), (24, 27)],
                   Feature(gene='sagD', location='13'): [(13, 15)]
               })
        interval_metadata1.concat(interval_metadata2, inplace=True)
        self.assertEqual(interval_metadata1.features, {
                         Feature(gene='sagA', location='0'): [(0, 2),
                                                              (4, 7)],
                         Feature(gene='sagB', location='3'): [(3, 5)],
                         Feature(gene='sagC', location='10'): [(10, 12),
                                                               (24, 27)],
                         Feature(gene='sagD', location='13'): [(13, 15)]})

    def test_eq(self):
        interval_metadata1 = IntervalMetadata(features={
                   Feature(gene='sagA', location='0'): [(0, 2), (4, 7)],
                   Feature(gene='sagB', location='3'): [(3, 5)]
               })

        interval_metadata2 = IntervalMetadata(features={
                   Feature(gene='sagA', location='0'): [(0, 2), (4, 7)],
                   Feature(gene='sagB', location='3'): [(3, 5)]
               })
        self.assertTrue(interval_metadata1 == interval_metadata2)


if __name__ == '__main__':
    unittest.main()
