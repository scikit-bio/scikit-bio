# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from copy import deepcopy

from skbio.metadata._interval import (_assert_valid_location,
                                      _assert_valid_boundary)
from skbio.metadata import Interval
from skbio.metadata import IntervalMetadata
from skbio.metadata._intersection import IntervalTree


class TestInterval(unittest.TestCase):
    def test_init_default(self):
        f = Interval(interval_metadata=IntervalMetadata(),
                     locations=[(1, 2), (4, 7)])

        self.assertTrue(f._interval_metadata is not None)
        self.assertListEqual(f.locations, [(1, 2), (4, 7)])
        self.assertListEqual(f.boundaries, [(True, True), (True, True)])
        self.assertDictEqual(f.metadata, {})

    def test_init(self):
        f = Interval(interval_metadata=IntervalMetadata(),
                     locations=[(1, 2), (4, 7)],
                     boundaries=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})

        self.assertTrue(f._interval_metadata is not None)
        self.assertListEqual(f.locations, [(1, 2), (4, 7)])
        self.assertListEqual(f.boundaries, [(True, False), (False, False)])
        self.assertDictEqual(f.metadata, {'name': 'sagA',
                                          'function': 'transport'})

    def test_init_iterables(self):
        f = Interval(interval_metadata=IntervalMetadata(),
                     locations=((1, 2), (4, 7)),
                     boundaries=((True, False), (False, False)),
                     metadata={'name': 'sagA', 'function': 'transport'})

        self.assertTrue(f._interval_metadata is not None)
        self.assertListEqual(f.locations, [(1, 2), (4, 7)])
        self.assertListEqual(f.boundaries, [(True, False), (False, False)])
        self.assertDictEqual(f.metadata, {'name': 'sagA',
                                          'function': 'transport'})

    def test_init_generator(self):
        def gen():
            for x in [(1, 2), (4, 7)]:
                yield x

        f = Interval(interval_metadata=IntervalMetadata(),
                     locations=gen(),
                     boundaries=((True, False), (False, False)),
                     metadata={'name': 'sagA', 'function': 'transport'})

        self.assertTrue(f._interval_metadata is not None)
        self.assertListEqual(f.locations, [(1, 2), (4, 7)])
        self.assertListEqual(f.boundaries, [(True, False), (False, False)])
        self.assertDictEqual(f.metadata, {'name': 'sagA',
                                          'function': 'transport'})

    def test_init_locations_scrambled(self):
        f = Interval(interval_metadata=IntervalMetadata(),
                     locations=[(4, 7), (1, 2)],
                     boundaries=[(True, False), (False, True)],
                     metadata={'name': 'sagA', 'function': 'transport'})

        self.assertTrue(f._interval_metadata is not None)
        self.assertListEqual(f.locations, [(1, 2), (4, 7)])
        self.assertListEqual(f.boundaries, [(False, True), (True, False)])
        self.assertDictEqual(f.metadata, {'name': 'sagA',
                                          'function': 'transport'})

    def test_init_bad(self):
        with self.assertRaises(TypeError):
            Interval(interval_metadata=None,
                     locations=[(4, 7)],
                     metadata={'name': 'sagA', 'function': 'transport'})

    def test_init_bad_locations(self):
        with self.assertRaises(TypeError):
            Interval(interval_metadata=IntervalMetadata(),
                     locations=[1, (4, 7)],
                     boundaries=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})

    def test_init_bad_boundaries(self):
        with self.assertRaises(ValueError):
            Interval(interval_metadata=IntervalMetadata(),
                     locations=[(1, 2), (4, 7)],
                     boundaries=[(False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})

    def test_repr(self):
        f = Interval(interval_metadata=IntervalMetadata(),
                     locations=[(1, 2)],
                     metadata={'name': 'sagA'})
        exp = (r"Interval\(interval_metadata=<[0-9]+>, locations=\[\(1, 2\)\],"
               " boundaries=\[\(True, True\)\], metadata={'name': 'sagA'}\)")
        obs = repr(f)
        self.assertRegex(obs, exp)
        # test for dropped
        f.drop()
        exp = (r"Interval\(dropped=True, locations=\[\(1, 2\)\],"
               " boundaries=\[\(True, True\)\], metadata={'name': 'sagA'}\)")
        obs = repr(f)
        self.assertRegex(obs, exp)

    def test_drop(self):
        f = Interval(interval_metadata=IntervalMetadata(),
                     locations=[(1, 2)],
                     metadata={'name': 'sagA'})
        f.drop()
        self.assertTrue(f._interval_metadata is None)
        self.assertTrue(f.dropped)
        self.assertTrue(f.locations, [(1, 2)])
        self.assertTrue(f.metadata, {'name': 'sagA'})
        # test the idempotence
        f.drop()
        self.assertTrue(f._interval_metadata is None)
        self.assertTrue(f.dropped)
        self.assertTrue(f.locations, [(1, 2)])
        self.assertTrue(f.metadata, {'name': 'sagA'})

    def test_equal(self):
        f1 = Interval(interval_metadata=IntervalMetadata(),
                      locations=[(1, 2), (4, 7)],
                      boundaries=[(True, False), (False, False)],
                      metadata={'name': 'sagA', 'function': 'transport'})

        f2 = Interval(interval_metadata=IntervalMetadata(),
                      locations=[(1, 2), (4, 7)],
                      boundaries=[(True, False), (False, False)],
                      metadata={'name': 'sagA', 'function': 'transport'})

        f3 = Interval(interval_metadata=IntervalMetadata(),
                      locations=[(1, 2), (4, 8)],
                      boundaries=[(True, True), (False, False)],
                      metadata={'name': 'sagA', 'function': 'transport'})

        f4 = Interval(interval_metadata=IntervalMetadata(),
                      locations=[(1, 2), (4, 8)],
                      boundaries=[(True, False), (False, False)],
                      metadata={'name': 'sagB', 'function': 'transport'})
        self.assertEqual(f2, f1)
        self.assertNotEqual(f1, f3)
        self.assertNotEqual(f1, f4)
        self.assertNotEqual(f4, f3)

    def test_equal_scrambled(self):
        im = IntervalMetadata()
        f1 = Interval(locations=[(9, 12), (4, 5)],
                      metadata={'name': 'sagA', 'function': 'transport'},
                      interval_metadata=im)
        f2 = Interval(locations=[(4, 5), (9, 12)],
                      metadata={'name': 'sagA', 'function': 'transport'},
                      interval_metadata=im)
        self.assertEqual(f1, f2)

    def test_get_locations(self):
        im = IntervalMetadata()
        f = Interval(interval_metadata=im,
                     locations=[(1, 2), (4, 7)],
                     boundaries=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        self.assertEqual(f.locations, [(1, 2), (4, 7)])
        self.assertEqual(im._is_stale_tree, True)

    def test_set_locations(self):
        im = IntervalMetadata()
        f = Interval(interval_metadata=im,
                     locations=[(1, 2), (4, 7)],
                     boundaries=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        f.locations = [(4, 7), (1, 3)]
        self.assertEqual(f.locations, [(1, 3), (4, 7)])
        self.assertEqual(f.boundaries, [(True, True), (True, True)])
        self.assertEqual(im._is_stale_tree, True)

    def test_set_locations_bad(self):
        f = Interval(interval_metadata=IntervalMetadata(),
                     locations=[(1, 2), (4, 7)],
                     boundaries=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        for value in [1, 's']:
            with self.assertRaises(TypeError):
                f.locations = value

        for value in [[(3, 1)], [('s', 1)], ()]:
            with self.assertRaises(ValueError):
                f.locations = value

    def test_get_boundaries(self):
        im = IntervalMetadata()
        f = Interval(interval_metadata=im,
                     locations=[(1, 2), (4, 7)],
                     boundaries=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        self.assertEqual(f.boundaries, [(True, False), (False, False)])

    def test_set_boundaries(self):
        im = IntervalMetadata()
        f = Interval(interval_metadata=im,
                     locations=[(1, 2), (4, 7)],
                     boundaries=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        f.boundaries = [(False, False), (False, False)]
        self.assertEqual(f.boundaries, [(False, False), (False, False)])

    def test_set_boundaries_bad(self):
        im = IntervalMetadata()
        f = Interval(interval_metadata=im,
                     locations=[(1, 2), (4, 7)],
                     boundaries=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        for value in [[(False, False)], ()]:
            with self.assertRaises(ValueError):
                f.boundaries = value
        for value in [1, True]:
            with self.assertRaises(TypeError):
                f.boundaries = value

    def test_get_metadata(self):
        im = IntervalMetadata()
        f = Interval(interval_metadata=im,
                     locations=[(1, 2), (4, 7)],
                     boundaries=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        f.metadata['name'] = 'sagB'
        self.assertEqual(f.metadata, {'name': 'sagB', 'function': 'transport'})

    def test_set_metadata(self):
        im = IntervalMetadata()
        f = Interval(interval_metadata=im,
                     locations=[(1, 2), (4, 7)],
                     boundaries=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        f.metadata = {'name': 'sagB', 'function': 'transport'}
        self.assertDictEqual(f.metadata,
                             {'name': 'sagB', 'function': 'transport'})
        f.metadata = {}
        self.assertDictEqual(f.metadata, {})

    def test_set_metadata_bad(self):
        im = IntervalMetadata()
        f = Interval(interval_metadata=im,
                     locations=[(1, 2), (4, 7)],
                     boundaries=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        for value in [1, '', None]:
            with self.assertRaises(TypeError):
                f.metadata = value

    def test_set_on_dropped(self):
        im = IntervalMetadata()
        f = Interval(interval_metadata=im,
                     locations=[(1, 2)],
                     boundaries=[(True, False)],
                     metadata={'name': 'sagA'})
        f.drop()
        with self.assertRaises(RuntimeError):
            f.boundaries = None
        with self.assertRaises(RuntimeError):
            f.locations = [(1, 2)]
        with self.assertRaises(RuntimeError):
            f.metadata = {}

    def test_get_on_dropped(self):
        im = IntervalMetadata()
        f = Interval(interval_metadata=im,
                     locations=[(1, 2)],
                     boundaries=[(True, False)],
                     metadata={'name': 'sagA'})
        f.drop()

        self.assertEqual(f.boundaries, [(True, False)])
        self.assertEqual(f.locations, [(1, 2)])
        self.assertEqual(f.metadata, {'name': 'sagA'})


class TestIntervalUtil(unittest.TestCase):
    def test_assert_valid_location(self):
        intvls = [(1, 2), (-1, 2)]
        for intvl in intvls:
            try:
                _assert_valid_location(intvl)
            except TypeError:
                self.assertTrue(False)

    def test_assert_valid_location_wrong_type(self):
        intvls = [[1, 2], 1, [1, 2, 3]]
        for intvl in intvls:
            with self.assertRaises(TypeError):
                _assert_valid_location(intvl)

    def test_assert_valid_location_wrong_value(self):
        intvls = [(1, 2, 3), (2, 1), (True, 0), ('s', 'r')]
        for intvl in intvls:
            with self.assertRaises(ValueError):
                _assert_valid_location(intvl)

    def test_assert_valid_boundary(self):
        boundaries = [(True, False), (True, True)]
        for boundary in boundaries:
            try:
                _assert_valid_boundary(boundary)
            except:
                self.assertTrue(False)

    def test_assert_valid_boundary_wrong_value(self):
        boundaries = [(True, False, True), ()]
        for boundary in boundaries:
            with self.assertRaises(ValueError):
                _assert_valid_boundary(boundary)

    def test_assert_valid_boundary_wrong_type(self):
        boundaries = [[True, False], 's', 1, (0, 1), ('s', '')]
        for boundary in boundaries:
            with self.assertRaises(TypeError):
                _assert_valid_boundary(boundary)


class TestIntervalMetadata(unittest.TestCase):
    def setUp(self):
        self.im_empty = IntervalMetadata()
        self.im_1 = IntervalMetadata()
        self.interval_1 = Interval(
            interval_metadata=self.im_1,
            locations=[(1, 2), (4, 7)],
            metadata={'gene': 'sagA',  'location': 0})
        self.im_2 = deepcopy(self.im_1)
        self.interval_2 = Interval(
            interval_metadata=self.im_2,
            locations=[(3, 5)],
            metadata={'gene': 'sagB', 'location': 0})

    def test_init(self):
        self.assertFalse(self.im_empty._is_stale_tree)
        self.assertEqual(self.im_empty._intervals, [])

    def test_duplicate(self):
        '''Test query and drop methods on duplicate Intervals.'''
        intvl_1 = self.im_empty.add([(1, 2)])
        intvl_2 = self.im_empty.add([(1, 2)])
        self.assertEqual(len(list(self.im_empty.query([(1, 2)]))), 2)
        self.im_empty.drop([intvl_1])
        self.assertEqual(len(self.im_empty._intervals), 1)
        self.assertTrue(self.im_empty._intervals[0] is intvl_2)

    def test_duplicate_locations(self):
        intvl = self.im_empty.add([(1, 2), (1, 2)])
        intvls = list(self.im_empty.query([(1, 2)]))
        self.assertEqual(len(intvls), 1)
        self.assertTrue(intvl is intvls[0])

    def test_sort(self):
        interval = Interval(
            self.im_2,
            [(1, 2), (3, 8)],
            metadata={'gene': 'sagA',  'location': 0})
        im = deepcopy(self.im_2)
        self.im_2.sort(False)
        # check sorting does not have other side effects
        self.assertEqual(im, self.im_2)
        self.assertEqual(self.im_2._intervals,
                         [self.interval_2, interval, self.interval_1])

        self.im_2.sort()
        self.assertEqual(im, self.im_2)
        self.assertEqual(self.im_2._intervals,
                         [self.interval_1, interval, self.interval_2])

        self.im_empty.sort()
        self.assertEqual(self.im_empty, IntervalMetadata())

    def test_add(self):
        self.im_empty.add(locations=[(1, 2), (4, 7)],
                          metadata={'gene': 'sagA',  'location': 0})
        self.assertTrue(self.im_empty._is_stale_tree)
        interval = self.im_empty._intervals[0]
        self.assertEqual(interval.locations, [(1, 2), (4, 7)])
        self.assertEqual(interval.metadata, {'gene': 'sagA', 'location': 0})
        self.assertTrue(isinstance(self.im_empty._interval_tree, IntervalTree))

    def test_query_attribute(self):
        intervals = self.im_2._query_attribute({})
        for i, j in zip(intervals, self.im_2._intervals):
            self.assertEqual(i, j)

        intervals = list(self.im_2._query_attribute(None))
        self.assertEqual(len(intervals), 0)

        for i in self.im_2._intervals:
            intervals = list(self.im_2._query_attribute(i.metadata))
            self.assertEqual(len(intervals), 1)
            self.assertEqual(intervals[0], i)

    def test_query_interval(self):
        intervals = list(self.im_2._query_interval((1, 2)))
        self.assertEqual(len(intervals), 1)
        self.assertEqual(intervals[0], self.interval_1)

        intervals = list(self.im_2._query_interval((3, 4)))
        self.assertEqual(len(intervals), 1)
        self.assertEqual(intervals[0], self.interval_2)

        intervals = {repr(i) for i in self.im_2._query_interval((1, 7))}
        self.assertEqual(len(intervals), 2)
        self.assertSetEqual(intervals,
                            {repr(i) for i in self.im_2._intervals})

    def test_query(self):
        intervals = list(self.im_2.query(locations=[(1, 5)],
                                         metadata={'gene': 'sagA'}))
        self.assertEqual(len(intervals), 1)
        self.assertEqual(intervals[0], self.interval_1)

    def test_query_empty(self):
        intervals = list(self.im_1.query())
        self.assertEqual(len(intervals), 0)

    def test_query_no_hits(self):
        intervals = list(self.im_2.query(locations=[(100, 200)]))
        self.assertEqual(len(intervals), 0)

        intervals = list(self.im_2.query(metadata={'gene': 'sagC'}))
        self.assertEqual(len(intervals), 0)

        intervals = list(self.im_2.query(locations=[(1, 2)],
                                         metadata={'gene': 'sagC'}))
        self.assertEqual(len(intervals), 0)

    def test_query_interval_only(self):
        for loc in [[(1, 7)],
                    [(1, 2), (3, 4)]]:
            intervals = list(self.im_2.query(locations=loc))
            self.assertEqual(len(intervals), 2)
            self.assertEqual(intervals[0], self.interval_1)
            self.assertEqual(intervals[1], self.interval_2)

    def test_query_metadata_only(self):
        intervals = list(self.im_2.query(metadata={'gene': 'sagB'}))
        self.assertEqual(len(intervals), 1)
        self.assertEqual(intervals[0], self.interval_2)

        intervals = list(self.im_2.query(metadata={'location': 0}))
        self.assertEqual(len(intervals), 2)
        self.assertEqual(intervals[0], self.interval_1)
        self.assertEqual(intervals[1], self.interval_2)

    def test_query_negative(self):
        intervals = list(self.im_2.query(locations=[(100, 101)]))
        self.assertEqual(len(intervals), 0)

    def test_drop(self):
        intvl = self.im_2._intervals[0]
        self.im_2.drop([intvl])
        self.assertEqual(len(self.im_2._intervals), 1)
        self.assertEqual(self.im_2._intervals[0], self.interval_2)
        # test the intvl was set to dropped
        self.assertTrue(intvl.dropped)

    def test_drop_all(self):
        self.im_2.drop(self.im_2._intervals)
        self.assertEqual(self.im_2, self.im_empty)

    def test_reverse(self):
        self.im_2._reverse(length=10)
        Interval(
            interval_metadata=self.im_empty,
            locations=[(3, 6), (8, 9)],
            metadata={'gene': 'sagA',  'location': 0})
        Interval(
            interval_metadata=self.im_empty,
            locations=[(5, 7)],
            metadata={'gene': 'sagB', 'location': 0})
        self.assertEqual(self.im_2, self.im_empty)

    def test_eq(self):
        im1 = IntervalMetadata()
        im1.add(metadata={'gene': 'sagA', 'location': '0'},
                locations=[(0, 2), (4, 7)])
        im1.add(metadata={'gene': 'sagB', 'location': '3'},
                locations=[(3, 5)])

        # The ordering shouldn't matter
        im2 = IntervalMetadata()
        im2.add(metadata={'gene': 'sagB', 'location': '3'},
                locations=[(3, 5)])
        im2.add(metadata={'gene': 'sagA', 'location': '0'},
                locations=[(0, 2), (4, 7)])

        im3 = IntervalMetadata()
        im3.add(metadata={'gene': 'sagA', 'location': '3'},
                locations=[(0, 2), (4, 7)])
        im3.add(metadata={'gene': 'sagB', 'location': '3'},
                locations=[(3, 5)])

        self.assertEqual(im1, im2)
        self.assertNotEqual(im1, im3)

    def test_repr(self):
        exp = '''0 interval features
-------------------'''
        self.assertEqual(repr(self.im_empty), exp)

        self.im_empty.add([(1, 2)], metadata={'gene': 'sagA'})
        self.im_empty.add([(3, 4)], metadata={'gene': 'sagB'})
        exp = '''2 interval features
-------------------
Interval\(interval_metadata=<[0-9]+>, locations=\[\(1, 2\)\], \
boundaries=\[\(True, True\)\], metadata={'gene': 'sagA'}\)
Interval\(interval_metadata=<[0-9]+>, locations=\[\(3, 4\)\], \
boundaries=\[\(True, True\)\], metadata={'gene': 'sagB'}\)'''
        self.assertRegex(repr(self.im_empty), exp)

        self.im_empty.add([(3, 4)], metadata={'gene': 'sagC'})
        self.im_empty.add([(3, 4)], metadata={'gene': 'sagD'})
        self.im_empty.add([(3, 4)], metadata={'gene': 'sagE'})
        self.im_empty.add([(3, 4)], metadata={'gene': 'sagF'})
        exp = '''6 interval features
-------------------
Interval\(interval_metadata=<[0-9]+>, locations=\[\(1, 2\)\], \
boundaries=\[\(True, True\)\], metadata={'gene': 'sagA'}\)
Interval\(interval_metadata=<[0-9]+>, locations=\[\(3, 4\)\], \
boundaries=\[\(True, True\)\], metadata={'gene': 'sagB'}\)
...
Interval\(interval_metadata=<[0-9]+>, locations=\[\(3, 4\)\], \
boundaries=\[\(True, True\)\], metadata={'gene': 'sagE'}\)
Interval\(interval_metadata=<[0-9]+>, locations=\[\(3, 4\)\], \
boundaries=\[\(True, True\)\], metadata={'gene': 'sagF'}\)'''
        self.assertRegex(repr(self.im_empty), exp)


if __name__ == '__main__':
    unittest.main()
