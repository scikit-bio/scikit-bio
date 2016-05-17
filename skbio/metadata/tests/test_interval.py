# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import unittest

from skbio.metadata._interval import _assert_valid_interval
from skbio.metadata import Interval
from skbio.metadata import IntervalMetadata


class TestInterval(unittest.TestCase):

    def test_constructor(self):
        f = Interval(interval_metadata=IntervalMetadata(),
                     intervals=[(1, 2), (4, 7)],
                     boundaries=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})

        self.assertTrue(f._interval_metadata is not None)
        self.assertListEqual(f.intervals, [(1, 2), (4, 7)])
        self.assertListEqual(f.boundaries, [(True, False), (False, False)])
        self.assertDictEqual(f.metadata, {'name': 'sagA',
                                          'function': 'transport'})

    def test_bad_constructor(self):
        with self.assertRaises(TypeError):
            Interval(interval_metadata=IntervalMetadata(),
                     intervals=[1, (4, 7)],
                     boundaries=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})

    def test_repr(self):
        im = IntervalMetadata()
        f = Interval(interval_metadata=im,
                     intervals=[(1, 2), (4, 7)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        exp = (r"<Interval: parent=%s, 2 intervals,"
               " 2 fields, dropped=False>") % (hex(id(im)))
        res = repr(f)
        self.assertEqual(res, exp)

    def test_cmp(self):
        f1 = Interval(interval_metadata=IntervalMetadata(),
                      intervals=[(1, 2)],
                      boundaries=[(True, False), (False, False)],
                      metadata={'name': 'sagA', 'function': 'transport'})
        f2 = Interval(interval_metadata=IntervalMetadata(),
                      intervals=[(10, 20)],
                      boundaries=[(True, False), (False, False)],
                      metadata={'name': 'sagA', 'function': 'transport'})
        self.assertTrue(f1._cmp(f2))
        self.assertFalse(f2._cmp(f1))

    def test_contains(self):
        f1 = Interval(interval_metadata=IntervalMetadata(),
                      intervals=[(1, 2)],
                      boundaries=[(True, False), (False, False)],
                      metadata={'name': 'sagA', 'function': 'transport'})
        self.assertIn('name', f1.metadata)
        self.assertIn('function', f1.metadata)
        self.assertFalse('gene' in f1.metadata)

    def test_hash(self):

        f1 = Interval(interval_metadata=IntervalMetadata(),
                      intervals=[(1, 2)],
                      boundaries=[(True, False), (False, False)],
                      metadata={'name': 'sagA', 'function': 'transport'})
        f2 = Interval(interval_metadata=IntervalMetadata(),
                      intervals=[(1, 2)],
                      boundaries=[(True, False), (False, False)],
                      metadata={'name': 'sagA', 'function': 'transport'})
        self.assertEqual(f1._hash(), f2._hash())

    def test_drop(self):
        im = IntervalMetadata()
        im.add(metadata={'gene': 'sagA', 'location': 0},
               intervals=[(0, 2), (4, 7)])
        self.assertEqual(im._metadata[0].intervals,
                         [(0, 2), (4, 7)])
        self.assertEqual(im._metadata[0].metadata,
                         {'gene': 'sagA', 'location': 0})
        self.assertTrue(im._intervals is not None)

        feats = im.query(intervals=[(1, 2)])

        self.assertEqual(len(feats), 1)
        self.assertEqual(feats[0].metadata, {'gene': 'sagA', 'location': 0})

        feats[0].drop()
        feats = im.query(intervals=[(1, 2)])
        self.assertEqual(len(feats), 0)

    def test_drop_none(self):
        im = IntervalMetadata()
        im.add(metadata={'gene': 'sagA', 'location': 0},
               intervals=[(0, 2), (4, 7)])
        invs = im.query(intervals=[(0, 7)])
        im.drop(intervals=[(0, 7)])
        invs[0].drop()

        feats = im.query(intervals=[(1, 2)])
        self.assertEqual(len(feats), 0)

    def test_equal(self):
        f = Interval(interval_metadata=IntervalMetadata(),
                     intervals=[(1, 2), (4, 7)],
                     boundaries=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})

        f1 = Interval(interval_metadata=IntervalMetadata(),
                      intervals=[(1, 2), (4, 7)],
                      boundaries=[(True, False), (False, False)],
                      metadata={'name': 'sagA', 'function': 'transport'})

        f2 = Interval(interval_metadata=IntervalMetadata(),
                      intervals=[(1, 2), (4, 8)],
                      boundaries=[(True, False), (False, False)],
                      metadata={'name': 'sagA', 'function': 'transport'})

        f3 = Interval(interval_metadata=IntervalMetadata(),
                      intervals=[(1, 2), (4, 8)],
                      boundaries=[(True, False), (False, False)],
                      metadata={'name': 'sagB', 'function': 'transport'})
        self.assertEqual(f, f1)
        self.assertNotEqual(f, f2)
        self.assertNotEqual(f, f3)
        self.assertNotEqual(f2, f3)

    def test_equal_scrambled(self):
        im = IntervalMetadata()
        f1 = Interval(intervals=[(9, 12), (4, 5)],
                      metadata={'name': 'sagA', 'function': 'transport'},
                      interval_metadata=im)
        f2 = Interval(intervals=[(4, 5), (9, 12)],
                      metadata={'name': 'sagA', 'function': 'transport'},
                      interval_metadata=im)
        self.assertEqual(f1, f2)

    def test_set_interval(self):
        im = IntervalMetadata()
        f = Interval(interval_metadata=im,
                     intervals=[(1, 2), (4, 7)],
                     boundaries=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        f.intervals = [(1, 3), (4, 7)]
        self.assertEqual(f.intervals, [(1, 3), (4, 7)])
        self.assertEqual(im._is_stale_tree, True)

    def test_set_interval_none(self):
        with self.assertRaises(RuntimeError):
            f = Interval(interval_metadata=None,
                         intervals=[(1, 2), (4, 7)],
                         boundaries=[(True, False), (False, False)],
                         metadata={'name': 'sagA', 'function': 'transport'})
            f.intervals = [(1, 3), (4, 7)]


class TestIntervalMetadata(unittest.TestCase):

    def test_init(self):
        im = IntervalMetadata()
        self.assertEqual(im._is_stale_tree, False)

    def test_assert_valid_interval_tuple(self):
        interval = (1, 2)
        _assert_valid_interval(interval)
        st, end = interval
        self.assertEqual(st, 1)
        self.assertEqual(end, 2)

    def test_assert_valid_interval_tuple_swapped(self):
        interval = (2, 1)
        with self.assertRaises(ValueError):
            _assert_valid_interval(interval)

    def test_assert_valid_interval_tuple_bad(self):
        with self.assertRaises(ValueError):
            _assert_valid_interval((1, 2, 3))

    def test_add(self):
        im = IntervalMetadata()
        im.add(intervals=[(1, 2), (4, 7)],
               metadata={'gene': 'sagA',  'location': 0})

        self.assertEqual(im._metadata[0].intervals,
                         [(1, 2), (4, 7)])
        self.assertEqual(im._metadata[0].metadata,
                         {'gene': 'sagA', 'location': 0})
        self.assertTrue(im._intervals is not None)

    def test_query(self):
        im = IntervalMetadata()
        im.add(intervals=[(0, 2), (4, 7)],
               metadata={'gene': 'sagA', 'location': 0})
        im.add(intervals=[(3, 5)],
               metadata={'gene': 'sagB', 'location': 0})

        feats = im.query(intervals=[(1, 2)])
        self.assertEqual(len(feats), 1)
        self.assertEqual(feats[0].metadata, {'gene': 'sagA', 'location': 0})

        feats = im.query(metadata={'gene': 'sagB'})
        self.assertEqual(len(feats), 1)
        self.assertEqual(feats[0].metadata, {'gene': 'sagB', 'location': 0})
        self.assertEqual(feats[0].intervals, [(3, 5)])

    def test_inconsistent_query(self):
        im = IntervalMetadata()
        im.add(intervals=[(3, 5)],
               metadata={'gene': 'sagA', 'software': 'qiime'})
        im.add(intervals=[(3, 5)],
               metadata={'gene': 'sagA', 'software': 'chiime'})
        im.add(intervals=[(3, 5)],
               metadata={'gene': 'sagA', 'software': 'keimei'})
        feats = im.query(metadata={'gene': 'sagA'})
        self.assertEqual(len(feats), 3)

        feats = im.query(metadata={'software': 'qiime'})
        self.assertEqual(len(feats), 1)
        self.assertEqual(feats[0].metadata,
                         {'gene': 'sagA', 'software': 'qiime'})
        self.assertEqual(feats[0].intervals, [(3, 5)])

    def test_query_duplicate1(self):
        im = IntervalMetadata()
        im.add(metadata={'gene': 'sagA', 'location': 0},
               intervals=[(0, 2), (4, 7)])
        im.add(metadata={'gene': 'sagB', 'location': 0},
               intervals=[(3, 5)])

        feats = im.query(intervals=[(1, 5)])
        exp = [Interval(metadata={'gene': 'sagA', 'location': 0},
                        intervals=[(0, 2), (4, 7)],
                        interval_metadata=im),
               Interval(metadata={'gene': 'sagB', 'location': 0},
                        intervals=[(3, 5)],
                        interval_metadata=im)]
        self.assertListEqual(sorted(feats, key=lambda x: x.intervals),
                             sorted(exp, key=lambda x: x.intervals))

    def test_query_negative(self):
        interval_metadata = IntervalMetadata()
        interval_metadata.add(intervals=[(40, 70)],
                              boundaries=None,
                              metadata={'name': 'sagB'})
        res = interval_metadata.query(intervals=[(100, 101)])
        self.assertEqual(len(res), 0)

    def test_query_intersection(self):
        im = IntervalMetadata()
        im.add(metadata={'gene': 'sagA', 'location': 0},
               intervals=[(0, 2), (4, 7)])
        im.add(metadata={'gene': 'sagA', 'location': 0},
               intervals=[(10, 15)])
        im.add(metadata={'gene': 'sagB', 'location': 0},
               intervals=[(3, 5)])

        feats = im.query(intervals=[(1, 5)], metadata={'gene': 'sagA'})
        exp = [Interval(metadata={'gene': 'sagA', 'location': 0},
                        intervals=[(0, 2), (4, 7)],
                        interval_metadata=im)]

        self.assertEqual(feats, exp)

    def test_set_interval_interval(self):
        interval_metadata = IntervalMetadata()
        interval_metadata.add(intervals=[(0, 2), (4, 7)],
                              boundaries=None, metadata={'name': 'sagA'})
        interval_metadata.add(intervals=[(40, 70)],
                              boundaries=None, metadata={'name': 'sagA'})
        interval_metadata.add(intervals=[(3, 4)],
                              boundaries=None, metadata={'name': 'sagB'})
        feats = list(interval_metadata.query(intervals=[(1, 2)]))
        feats[0].intervals = [(1, 2)]
        feats = list(interval_metadata.query(intervals=[(1, 2)]))
        self.assertEqual(feats[0].intervals, [(1, 2)])

    def test_set_interval_attribute(self):
        interval_metadata = IntervalMetadata()
        interval_metadata.add(intervals=[(0, 2), (4, 7)],
                              boundaries=None, metadata={'name': 'sagA'})
        interval_metadata.add(intervals=[(40, 70)],
                              boundaries=None, metadata={'name': 'sagA'})
        interval_metadata.add(intervals=[(3, 4)],
                              boundaries=None, metadata={'name': 'sagB'})
        feats = list(interval_metadata.query(intervals=[(1, 2)]))
        feats[0].metadata['name'] = 'sagC'
        feats = list(interval_metadata.query(intervals=[(1, 2)]))
        self.assertEqual(feats[0].metadata['name'], 'sagC')

    def test_drop(self):
        interval_metadata = IntervalMetadata()
        interval_metadata.add(intervals=[(0, 2), (4, 7)],
                              boundaries=None, metadata={'name': 'sagA'})
        interval_metadata.add(intervals=[(40, 70)],
                              boundaries=None, metadata={'name': 'sagA'})
        interval_metadata.add(intervals=[(3, 4)],
                              boundaries=None, metadata={'name': 'sagB'})
        interval_metadata.drop(metadata={'name': 'sagA'})
        feats = list(interval_metadata.query(metadata={'name': 'sagA'}))
        self.assertEqual(len(feats), 0)
        feats = list(interval_metadata.query(metadata={'name': 'sagB'}))
        self.assertGreater(len(feats), 0)

    def test_drop_interval(self):
        from skbio.metadata import IntervalMetadata
        interval_metadata = IntervalMetadata()
        interval_metadata.add(intervals=[(0, 2), (4, 7)],
                              boundaries=None, metadata={'name': 'sagA'})
        interval_metadata.add(intervals=[(40, 70)],
                              boundaries=None, metadata={'name': 'sagB'})
        interval_metadata.drop(metadata={'name': 'sagA'})
        res = interval_metadata.query(intervals=[(1, 2)])
        self.assertEqual(len(res), 0)
        feats = list(interval_metadata.query(metadata={'name': 'sagB'}))
        self.assertGreater(len(feats), 0)

    def test_reverse_complement(self):
        interval_metadata = IntervalMetadata()
        interval_metadata.add(metadata={'gene': 'sagB', 'location': 0},
                              intervals=[(3, 5)])
        interval_metadata._reverse(length=10)
        feats = interval_metadata.query([(5, 7)])
        exp = Interval(intervals=[(5, 7)],
                       metadata={'gene': 'sagB', 'location': 0},
                       interval_metadata=interval_metadata)
        self.assertEqual(feats[0], exp)

    def test_eq(self):
        interval_metadata1 = IntervalMetadata()
        interval_metadata1.add(metadata={'gene': 'sagA', 'location': '0'},
                               intervals=[(0, 2), (4, 7)])
        interval_metadata1.add(metadata={'gene': 'sagB', 'location': '3'},
                               intervals=[(3, 5)])

        interval_metadata2 = IntervalMetadata()
        interval_metadata2.add(metadata={'gene': 'sagA', 'location': '0'},
                               intervals=[(0, 2), (4, 7)])
        interval_metadata2.add(metadata={'gene': 'sagB', 'location': '3'},
                               intervals=[(3, 5)])

        interval_metadata3 = IntervalMetadata()
        interval_metadata3.add(metadata={'gene': 'sagA', 'location': '3'},
                               intervals=[(0, 2), (4, 7)])
        interval_metadata3.add(metadata={'gene': 'sagB', 'location': '3'},
                               intervals=[(3, 5)])

        # The ordering shouldn't matter
        interval_metadata4 = IntervalMetadata()
        interval_metadata4.add(metadata={'gene': 'sagB', 'location': '3'},
                               intervals=[(3, 5)])
        interval_metadata4.add(metadata={'gene': 'sagA', 'location': '0'},
                               intervals=[(0, 2), (4, 7)])

        self.assertEqual(interval_metadata1, interval_metadata2)
        self.assertNotEqual(interval_metadata1, interval_metadata3)
        self.assertEqual(interval_metadata1, interval_metadata4)

    def test_repr(self):

        interval_metadata = IntervalMetadata()
        interval_metadata.add(metadata={'gene': 'sagA', 'location': '0'},
                              intervals=[(0, 2), (4, 7)])
        interval_metadata.add(metadata={'gene': 'sagB', 'location': '3'},
                              intervals=[(3, 15)])
        interval_metadata.add(metadata={'gene': 'sagC', 'location': '3'},
                              intervals=[(4, 15)])
        interval_metadata.add(metadata={'gene': 'sagD', 'location': '3'},
                              intervals=[(5, 15)])
        interval_metadata.add(metadata={'gene': 'sagE', 'location': '3'},
                              intervals=[(6, 15)])
        interval_metadata.add(metadata={'gene': 'sagF', 'location': '3'},
                              intervals=[(7, 15)])
        interval_metadata.add(metadata={'gene': 'sagG', 'location': '3'},
                              intervals=[(8, 15)])
        interval_metadata.add(metadata={'gene': 'sagH', 'location': '3'},
                              intervals=[(9, 15)])
        interval_metadata.add(metadata={'gene': 'sagI', 'location': '3'},
                              intervals=[(10, 15)])
        interval_metadata.add(metadata={'gene': 'sagJ', 'location': '3'},
                              intervals=[(13, 15)])
        interval_metadata.add(metadata={'gene': 'sagK', 'location': '3'},
                              intervals=[(14, 15)])
        interval_metadata2 = IntervalMetadata()
        exp = [Interval(intervals=[(0, 2), (4, 7)],
                        metadata={'location': '0', 'gene': 'sagA'},
                        interval_metadata=interval_metadata2),
               Interval(intervals=[(3, 15)],
                        metadata={'location': '3', 'gene': 'sagB'},
                        interval_metadata=interval_metadata2),
               Interval(intervals=[(4, 15)],
                        metadata={'location': '3', 'gene': 'sagC'},
                        interval_metadata=interval_metadata2),
               Interval(intervals=[(5, 15)],
                        metadata={'location': '3', 'gene': 'sagD'},
                        interval_metadata=interval_metadata2),
               Interval(intervals=[(6, 15)],
                        metadata={'location': '3', 'gene': 'sagE'},
                        interval_metadata=interval_metadata2),
               '...',
               Interval(intervals=[(8, 15)],
                        metadata={'location': '3', 'gene': 'sagG'},
                        interval_metadata=interval_metadata2),
               Interval(intervals=[(9, 15)],
                        metadata={'location': '3', 'gene': 'sagH'},
                        interval_metadata=interval_metadata2),
               Interval(intervals=[(10, 15)],
                        metadata={'location': '3', 'gene': 'sagI'},
                        interval_metadata=interval_metadata2),
               Interval(intervals=[(13, 15)],
                        metadata={'location': '3', 'gene': 'sagJ'},
                        interval_metadata=interval_metadata2),
               Interval(intervals=[(14, 15)],
                        metadata={'location': '3', 'gene': 'sagK'},
                        interval_metadata=interval_metadata2)]
        self.assertEqual(str(exp), repr(interval_metadata))

if __name__ == '__main__':
    unittest.main()
