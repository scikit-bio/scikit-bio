# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from copy import deepcopy, copy

from skbio.metadata._interval import (_assert_valid_bound,
                                      _assert_valid_fuzzy)
from skbio.metadata import Interval
from skbio.metadata import IntervalMetadata
from skbio.metadata._intersection import IntervalTree
from skbio.util._testing import ReallyEqualMixin


class TestInterval(unittest.TestCase, ReallyEqualMixin):
    def setUp(self):
        self.upper_bound = 100
        self.im = IntervalMetadata(self.upper_bound)

    def test_init_default(self):
        f = Interval(self.im, bounds=[(0, 2), (4, self.upper_bound)])

        self.assertTrue(f._interval_metadata is not None)
        self.assertListEqual(f.bounds, [(0, 2), (4, self.upper_bound)])
        self.assertListEqual(f.fuzzy, [(False, False), (False, False)])
        self.assertDictEqual(f.metadata, {})

    def test_init(self):
        f = Interval(interval_metadata=self.im,
                     bounds=[(1, 2), (4, 7)],
                     fuzzy=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})

        self.assertTrue(f._interval_metadata is not None)
        self.assertListEqual(f.bounds, [(1, 2), (4, 7)])
        self.assertListEqual(f.fuzzy, [(True, False), (False, False)])
        self.assertDictEqual(f.metadata, {'name': 'sagA',
                                          'function': 'transport'})

    def test_init_iterables(self):
        f = Interval(interval_metadata=self.im,
                     bounds=((1, 2), (4, 7)),
                     fuzzy=((True, False), (False, False)),
                     metadata={'name': 'sagA', 'function': 'transport'})

        self.assertTrue(f._interval_metadata is not None)
        self.assertListEqual(f.bounds, [(1, 2), (4, 7)])
        self.assertListEqual(f.fuzzy, [(True, False), (False, False)])
        self.assertDictEqual(f.metadata, {'name': 'sagA',
                                          'function': 'transport'})

    def test_init_generator(self):
        def gen():
            for x in [(1, 2), (4, 7)]:
                yield x

        f = Interval(interval_metadata=self.im,
                     bounds=gen(),
                     fuzzy=((True, False), (False, False)),
                     metadata={'name': 'sagA', 'function': 'transport'})

        self.assertTrue(f._interval_metadata is not None)
        self.assertListEqual(f.bounds, [(1, 2), (4, 7)])
        self.assertListEqual(f.fuzzy, [(True, False), (False, False)])
        self.assertDictEqual(f.metadata, {'name': 'sagA',
                                          'function': 'transport'})

    def test_init_bounds_scrambled(self):
        f = Interval(interval_metadata=self.im,
                     bounds=[(4, 7), (1, 2)],
                     fuzzy=[(True, False), (False, True)],
                     metadata={'name': 'sagA', 'function': 'transport'})

        self.assertTrue(f._interval_metadata is not None)
        self.assertListEqual(f.bounds, [(1, 2), (4, 7)])
        self.assertListEqual(f.fuzzy, [(False, True), (True, False)])
        self.assertDictEqual(f.metadata, {'name': 'sagA',
                                          'function': 'transport'})

    def test_init_no_interval_metadata(self):
        with self.assertRaises(TypeError):
            Interval(interval_metadata=None,
                     bounds=[(4, 7)],
                     metadata={'name': 'sagA', 'function': 'transport'})

    def test_init_empty_metadata(self):
        for i in 0, 1:
            # test that no exception is raised
            Interval(interval_metadata=self.im, bounds=[(i, i)])

    def test_init_out_of_bounds(self):
        with self.assertRaises(ValueError):
            Interval(interval_metadata=self.im,
                     bounds=[(1, 2), (4, 101)],
                     fuzzy=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        with self.assertRaises(ValueError):
            Interval(interval_metadata=self.im,
                     bounds=[(-1, 2), (4, 6)],
                     fuzzy=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})

    def test_init_bad_bounds(self):
        with self.assertRaises(TypeError):
            Interval(interval_metadata=self.im,
                     bounds=[1, (4, 7)],
                     fuzzy=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})

    def test_init_bad_fuzzy(self):
        with self.assertRaises(ValueError):
            Interval(interval_metadata=self.im,
                     bounds=[(1, 2), (4, 7)],
                     fuzzy=[(False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})

    def test_repr(self):
        f = Interval(interval_metadata=self.im,
                     bounds=[(1, 2)],
                     metadata={'name': 'sagA'})
        exp = (r"Interval\(interval_metadata=<[0-9]+>, bounds=\[\(1, 2\)\],"
               r" fuzzy=\[\(False, False\)\], metadata={'name': 'sagA'}\)")
        obs = repr(f)
        self.assertRegex(obs, exp)
        # test for dropped
        f.drop()
        exp = (r"Interval\(dropped=True, bounds=\[\(1, 2\)\],"
               r" fuzzy=\[\(False, False\)\], metadata={'name': 'sagA'}\)")
        obs = repr(f)
        self.assertRegex(obs, exp)

    def test_drop(self):
        f = Interval(interval_metadata=self.im,
                     bounds=[(1, 2)],
                     metadata={'name': 'sagA'})
        f.drop()
        self.assertTrue(f._interval_metadata is None)
        self.assertTrue(f.dropped)
        self.assertTrue(f.bounds, [(1, 2)])
        self.assertTrue(f.metadata, {'name': 'sagA'})
        # test the idempotence
        f.drop()
        self.assertTrue(f._interval_metadata is None)
        self.assertTrue(f.dropped)
        self.assertTrue(f.bounds, [(1, 2)])
        self.assertTrue(f.metadata, {'name': 'sagA'})

    def test_eq(self):
        f0 = Interval(interval_metadata=self.im,
                      bounds=[(4, 7), (1, 2)],
                      fuzzy=[(False, False), (True, False)],
                      metadata={'name': 'sagA', 'function': 'transport'})

        f1 = Interval(interval_metadata=self.im,
                      bounds=[(1, 2), (4, 7)],
                      fuzzy=[(True, False), (False, False)],
                      metadata={'name': 'sagA', 'function': 'transport'})

        f2 = Interval(interval_metadata=self.im,
                      bounds=[(1, 2), (4, 7)],
                      fuzzy=[(True, False), (False, False)],
                      metadata={'name': 'sagA', 'function': 'transport'})

        f3 = Interval(interval_metadata=self.im,
                      bounds=[(1, 2), (4, 7)],
                      fuzzy=[(True, True), (False, False)],
                      metadata={'name': 'sagA', 'function': 'transport'})

        f4 = Interval(interval_metadata=self.im,
                      bounds=[(1, 2), (4, 8)],
                      fuzzy=[(True, False), (False, False)],
                      metadata={'name': 'sagA', 'function': 'transport'})

        f5 = Interval(interval_metadata=self.im,
                      bounds=[(1, 2), (4, 7)],
                      fuzzy=[(True, False), (False, False)],
                      metadata={'name': 'sagB', 'function': 'transport'})

        # scramble bounds/fuzzy
        self.assertReallyEqual(f0, f1)
        self.assertReallyEqual(f2, f1)
        # diff fuzzy
        self.assertReallyNotEqual(f1, f3)
        # diff bounds
        self.assertReallyNotEqual(f1, f4)
        # diff metadata
        self.assertReallyNotEqual(f1, f5)

    def test_get_bounds(self):
        f = Interval(interval_metadata=self.im,
                     bounds=[(1, 2), (4, 7)],
                     fuzzy=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        self.assertEqual(f.bounds, [(1, 2), (4, 7)])
        self.assertEqual(self.im._is_stale_tree, True)

    def test_set_bounds(self):
        f = Interval(interval_metadata=self.im,
                     bounds=[(1, 2), (4, 7)],
                     fuzzy=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        f.bounds = [(4, 7), (1, 3)]
        self.assertEqual(f.bounds, [(1, 3), (4, 7)])
        self.assertEqual(f.fuzzy, [(False, False), (False, False)])
        self.assertEqual(self.im._is_stale_tree, True)

    def test_set_bounds_bad(self):
        f = Interval(interval_metadata=self.im,
                     bounds=[(1, 2), (4, 7)],
                     fuzzy=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        for value in [1, 's']:
            with self.assertRaises(TypeError):
                f.bounds = value

        for value in [[(-1, 2)],   # start < lower_bound
                      [(1, 101)],  # end > upper_bound
                      [(3, 1)],    # start < end
                      [('s', 1)], (), None]:  # invalid values
            with self.assertRaises(ValueError):
                f.bounds = value

    def test_get_fuzzy(self):
        f = Interval(interval_metadata=self.im,
                     bounds=[(1, 2), (4, 7)],
                     fuzzy=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        self.assertEqual(f.fuzzy, [(True, False), (False, False)])

    def test_set_fuzzy(self):
        f = Interval(interval_metadata=self.im,
                     bounds=[(1, 2), (4, 7)],
                     fuzzy=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        f.fuzzy = [(False, False), (False, False)]
        self.assertEqual(f.fuzzy, [(False, False), (False, False)])

    def test_set_fuzzy_bad(self):
        f = Interval(interval_metadata=self.im,
                     bounds=[(1, 2), (4, 7)],
                     fuzzy=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        for value in [[(False, False)], (), None]:
            with self.assertRaises(ValueError):
                f.fuzzy = value
        for value in [1, True]:
            with self.assertRaises(TypeError):
                f.fuzzy = value

    def test_delete_fuzzy(self):
        f = Interval(interval_metadata=self.im,
                     bounds=[(1, 2), (4, 7)],
                     fuzzy=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        del f.fuzzy
        self.assertEqual(f.fuzzy, [(False, False), (False, False)])
        # delete again
        del f.fuzzy
        self.assertEqual(f.fuzzy, [(False, False), (False, False)])

    def test_get_metadata(self):
        f = Interval(interval_metadata=self.im,
                     bounds=[(1, 2), (4, 7)],
                     fuzzy=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        f.metadata['name'] = 'sagB'
        self.assertEqual(f.metadata, {'name': 'sagB', 'function': 'transport'})

    def test_set_metadata(self):
        f = Interval(interval_metadata=self.im,
                     bounds=[(1, 2), (4, 7)],
                     fuzzy=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        f.metadata = {'name': 'sagB', 'function': 'transport'}
        self.assertDictEqual(f.metadata,
                             {'name': 'sagB', 'function': 'transport'})
        f.metadata = {}
        self.assertDictEqual(f.metadata, {})

    def test_set_metadata_bad(self):
        f = Interval(interval_metadata=self.im,
                     bounds=[(1, 2), (4, 7)],
                     fuzzy=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        for value in [1, '', None]:
            with self.assertRaises(TypeError):
                f.metadata = value

    def test_delete_metadata(self):
        f = Interval(interval_metadata=self.im,
                     bounds=[(1, 2), (4, 7)],
                     fuzzy=[(True, False), (False, False)],
                     metadata={'name': 'sagA', 'function': 'transport'})
        del f.metadata
        self.assertEqual(f.metadata, {})

    def test_set_delete_on_dropped(self):
        f = Interval(interval_metadata=self.im,
                     bounds=[(1, 2)],
                     fuzzy=[(True, False)],
                     metadata={'name': 'sagA'})
        f.drop()
        with self.assertRaises(RuntimeError):
            f.fuzzy = None
        with self.assertRaises(RuntimeError):
            f.bounds = [(1, 2)]
        with self.assertRaises(RuntimeError):
            f.metadata = {}
        with self.assertRaises(RuntimeError):
            del f.fuzzy
        with self.assertRaises(RuntimeError):
            del f.metadata

    def test_get_on_dropped(self):
        f = Interval(interval_metadata=self.im,
                     bounds=[(1, 2)],
                     fuzzy=[(True, False)],
                     metadata={'name': 'sagA'})
        f.drop()

        self.assertEqual(f.fuzzy, [(True, False)])
        self.assertEqual(f.bounds, [(1, 2)])
        self.assertEqual(f.metadata, {'name': 'sagA'})


class TestIntervalUtil(unittest.TestCase):
    def test_assert_valid_bound(self):
        intvls = [(1, 2), (-1, 2)]
        for intvl in intvls:
            try:
                _assert_valid_bound(intvl)
            except TypeError:
                self.assertTrue(False)

    def test_assert_valid_bound_wrong_type(self):
        intvls = [[1, 2], 1, [1, 2, 3]]
        for intvl in intvls:
            with self.assertRaises(TypeError):
                _assert_valid_bound(intvl)

    def test_assert_valid_bound_wrong_value(self):
        intvls = [(1, 2, 3), (2, 1), (True, 0), ('s', 'r')]
        for intvl in intvls:
            with self.assertRaises(ValueError):
                _assert_valid_bound(intvl)

    def test_assert_valid_fuzzy(self):
        fuzzy = [(True, False), (True, True)]
        for fuzzy in fuzzy:
            try:
                _assert_valid_fuzzy(fuzzy)
            except Exception:
                self.assertTrue(False)

    def test_assert_valid_fuzzy_wrong_value(self):
        fuzzy = [(True, False, True), ()]
        for fuzzy in fuzzy:
            with self.assertRaises(ValueError):
                _assert_valid_fuzzy(fuzzy)

    def test_assert_valid_fuzzy_wrong_type(self):
        fuzzy = [[True, False], 's', 1, (0, 1), ('s', '')]
        for fuzzy in fuzzy:
            with self.assertRaises(TypeError):
                _assert_valid_fuzzy(fuzzy)


class TestIntervalMetadata(unittest.TestCase, ReallyEqualMixin):
    def setUp(self):
        self.upper_bound = 10
        self.im_empty = IntervalMetadata(self.upper_bound)
        self.im_1 = IntervalMetadata(self.upper_bound)
        self.im_1_1 = Interval(
            interval_metadata=self.im_1,
            bounds=[(1, 2), (4, self.upper_bound)],
            metadata={'gene': 'sagA',  'bound': 0})
        self.im_2 = IntervalMetadata(self.upper_bound)
        self.im_2_1 = Interval(
            interval_metadata=self.im_2,
            bounds=[(1, 2), (4, self.upper_bound)],
            metadata={'gene': 'sagA',  'bound': 0})
        self.im_2_2 = Interval(
            interval_metadata=self.im_2,
            bounds=[(3, 5)],
            metadata={'gene': 'sagB', 'bound': 0, 'spam': [0]})

    def test_copy_empty(self):
        obs = copy(self.im_empty)
        self.assertEqual(obs, self.im_empty)
        self.assertIsNot(obs._intervals, self.im_empty._intervals)
        self.assertIsNot(obs._interval_tree, self.im_empty._interval_tree)

    def test_copy(self):
        obs = copy(self.im_2)
        self.assertEqual(obs, self.im_2)
        self.assertIsNot(obs._intervals, self.im_2._intervals)
        self.assertIsNot(obs._interval_tree, self.im_2._interval_tree)

        for i in range(self.im_2.num_interval_features):
            i1, i2 = obs._intervals[i], self.im_2._intervals[i]
            self.assertIsNot(i1, i2)
            self.assertIsNot(i1.bounds, i2.bounds)
            self.assertIsNot(i1.fuzzy, i2.fuzzy)
            self.assertIsNot(i1._interval_metadata, i2._interval_metadata)
            self.assertIsNot(i1.metadata, i2.metadata)
            for k in i1.metadata:
                self.assertIs(i1.metadata[k], i2.metadata[k])

    def test_deepcopy(self):
        obs = deepcopy(self.im_2)
        self.assertEqual(obs, self.im_2)
        self.assertIsNot(obs._intervals, self.im_2._intervals)
        self.assertIsNot(obs._interval_tree, self.im_2._interval_tree)

        for i in range(self.im_2.num_interval_features):
            i1, i2 = obs._intervals[i], self.im_2._intervals[i]
            self.assertIsNot(i1, i2)
            self.assertIsNot(i1.bounds, i2.bounds)
            self.assertIsNot(i1.fuzzy, i2.fuzzy)
            self.assertIsNot(i1.metadata, i2.metadata)

        i2.metadata['spam'].append(1)
        self.assertEqual(i2.metadata,
                         {'gene': 'sagB', 'bound': 0, 'spam': [0, 1]})
        self.assertEqual(i1.metadata,
                         {'gene': 'sagB', 'bound': 0, 'spam': [0]})

    def test_deepcopy_memo_is_respected(self):
        memo = {}
        deepcopy(self.im_1, memo)
        self.assertGreater(len(memo), 2)

    def test_init(self):
        self.assertFalse(self.im_empty._is_stale_tree)
        self.assertEqual(self.im_empty._intervals, [])

    def test_init_upper_bound_lt_lower_bound(self):
        # test that no exception is raised
        IntervalMetadata(0)

        with self.assertRaises(ValueError):
            IntervalMetadata(-1)

    def test_upper_bound_is_none(self):
        im = IntervalMetadata(None)
        # should not raise error
        im.add([(0, 1000000000)])
        self.assertIsNone(im.upper_bound)
        with self.assertRaisesRegex(
                TypeError, r'upper bound is `None`'):
            im._reverse()
        with self.assertRaisesRegex(
                TypeError, r'upper bound is `None`'):
            IntervalMetadata.concat([self.im_1, im])

    def test_init_copy_from(self):
        for i in [None, 99, 999]:
            obs = IntervalMetadata(i, self.im_1)
            exp = IntervalMetadata(i)
            exp.add(bounds=[(1, 2), (4, self.upper_bound)],
                    metadata={'gene': 'sagA',  'bound': 0})
            self.assertEqual(obs, exp)

    def test_init_copy_from_empty(self):
        for i in [None, 0, 9, 99, 999]:
            obs = IntervalMetadata(i, self.im_empty)
            exp = IntervalMetadata(i)
            self.assertEqual(obs, exp)
            # test it is shallow copy
            self.assertIsNot(obs._intervals, self.im_empty._intervals)
            self.assertIsNot(obs._interval_tree, self.im_empty._interval_tree)

    def test_init_copy_from_shallow_copy(self):
        obs = IntervalMetadata(self.upper_bound, self.im_2)
        self.assertEqual(self.im_2, obs)
        # test it is shallow copy
        self.assertIsNot(obs._intervals, self.im_2._intervals)
        self.assertIsNot(obs._interval_tree, self.im_2._interval_tree)
        for i in range(self.im_2.num_interval_features):
            i1, i2 = obs._intervals[i], self.im_2._intervals[i]
            self.assertIsNot(i1, i2)
            self.assertIsNot(i1.bounds, i2.bounds)
            self.assertIsNot(i1.fuzzy, i2.fuzzy)
            self.assertIsNot(i1._interval_metadata, i2._interval_metadata)
            self.assertIsNot(i1.metadata, i2.metadata)
            for k in i1.metadata:
                self.assertIs(i1.metadata[k], i2.metadata[k])

    def test_init_copy_from_error(self):
        i = self.upper_bound - 1
        with self.assertRaisesRegex(
                ValueError, r'larger than upper bound \(%r\)' % i):
            IntervalMetadata(i, self.im_2)

    def test_num_interval_features(self):
        self.assertEqual(self.im_empty.num_interval_features, 0)
        self.assertEqual(self.im_1.num_interval_features, 1)
        self.assertEqual(self.im_2.num_interval_features, 2)

    def test_duplicate(self):
        '''Test query and drop methods on duplicate Intervals.'''
        intvl_1 = self.im_empty.add([(1, 2)])
        intvl_2 = self.im_empty.add([(1, 2)])
        self.assertEqual(len(list(self.im_empty.query([(1, 2)]))), 2)
        self.im_empty.drop([intvl_1])
        self.assertEqual(len(self.im_empty._intervals), 1)
        self.assertTrue(self.im_empty._intervals[0] is intvl_2)

    def test_duplicate_bounds(self):
        intvl = self.im_empty.add([(1, 2), (1, 2)])
        intvls = list(self.im_empty.query([(1, 2)]))
        self.assertEqual(len(intvls), 1)
        self.assertTrue(intvl is intvls[0])

    def test_concat_empty(self):
        for i in 0, 1, 2:
            obs = IntervalMetadata.concat([self.im_empty] * i)
            exp = IntervalMetadata(self.upper_bound * i)
            self.assertEqual(obs, exp)

        obs = IntervalMetadata.concat([])
        self.assertEqual(obs, IntervalMetadata(0))

    def test_concat(self):
        im1 = IntervalMetadata(3)
        im2 = IntervalMetadata(4)
        im3 = IntervalMetadata(5)
        im1.add([(0, 2)], [(True, True)])
        im2.add([(0, 3)], [(True, False)], {'gene': 'sagA'})
        im2.add([(2, 4)], metadata={'gene': 'sagB'})
        im3.add([(1, 5)], [(False, True)], {'gene': 'sagC'})
        obs = IntervalMetadata.concat([im1, im2, im3])

        exp = IntervalMetadata(12)
        exp.add(bounds=[(0, 2)], fuzzy=[(True, True)])
        exp.add(bounds=[(3, 6)], fuzzy=[(True, False)],
                metadata={'gene': 'sagA'})
        exp.add(bounds=[(5, 7)], metadata={'gene': 'sagB'})
        exp.add(bounds=[(8, 12)], fuzzy=[(False, True)],
                metadata={'gene': 'sagC'})
        self.assertEqual(obs, exp)

    def test_merge(self):
        # empty + empty
        im = IntervalMetadata(self.upper_bound)
        self.im_empty.merge(im)
        self.assertEqual(self.im_empty, im)
        # empty + non-empty
        self.im_empty.merge(self.im_1)
        self.assertEqual(self.im_empty, self.im_1)
        # non-empty + non-empty
        self.im_empty.merge(self.im_2)
        self.im_2.merge(self.im_1)
        self.assertEqual(self.im_empty, self.im_2)

    def test_merge_unequal_upper_bounds(self):
        n = 3
        im1 = IntervalMetadata(n)
        for im in [self.im_empty, self.im_1]:
            with self.assertRaisesRegex(
                    ValueError,
                    r'not equal \(%d != %d\)' % (self.upper_bound, n)):
                im.merge(im1)

    def test_merge_to_unbounded(self):
        for im in [self.im_empty, self.im_1, IntervalMetadata(None)]:
            obs = IntervalMetadata(None)
            obs.merge(im)
            self.assertIsNone(obs.upper_bound)
            self.assertEqual(obs._intervals, im._intervals)

    def test_merge_unbounded_to_bounded(self):
        im = IntervalMetadata(None)
        with self.assertRaisesRegex(
                ValueError,
                r'Cannot merge an unbound IntervalMetadata object '
                'to a bounded one'):
            self.im_1.merge(im)
        # original im is not changed
        self.assertIsNone(im.upper_bound)
        self.assertEqual(im._intervals, [])

    def test_sort(self):
        interval = Interval(
            self.im_2,
            [(1, 2), (3, 8)],
            metadata={'gene': 'sagA',  'bound': 0})
        im = deepcopy(self.im_2)
        self.im_2.sort(False)
        # check sorting does not have other side effects
        self.assertEqual(im, self.im_2)
        self.assertEqual(self.im_2._intervals,
                         [self.im_2_2, self.im_2_1, interval])

        self.im_2.sort()
        self.assertEqual(im, self.im_2)
        self.assertEqual(self.im_2._intervals,
                         [interval, self.im_2_1, self.im_2_2])

        self.im_empty.sort()
        self.assertEqual(self.im_empty, IntervalMetadata(self.upper_bound))

    def test_add_eq_upper_bound(self):
        self.im_empty.add(bounds=[(1, 2), (4, self.upper_bound)],
                          metadata={'gene': 'sagA',  'bound': 0})
        self.assertTrue(self.im_empty._is_stale_tree)
        interval = self.im_empty._intervals[0]
        self.assertEqual(interval.bounds, [(1, 2), (4, self.upper_bound)])
        self.assertEqual(interval.metadata, {'gene': 'sagA', 'bound': 0})
        self.assertTrue(isinstance(self.im_empty._interval_tree, IntervalTree))

    def test_add_gt_upper_bound(self):
        with self.assertRaises(ValueError):
            self.im_empty.add(bounds=[(1, 2), (4, self.upper_bound+1)],
                              metadata={'gene': 'sagA',  'bound': 0})

    def test_add_eq_start_end_bound(self):
        for i in 0, 1, self.upper_bound:
            # test that no exception is raised
            self.im_empty.add(bounds=[(i, i)],
                              metadata={'gene': 'sagA',  'bound': 0})

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
        self.assertEqual(intervals[0], self.im_2_1)

        intervals = list(self.im_2._query_interval((3, 4)))
        self.assertEqual(len(intervals), 1)
        self.assertEqual(intervals[0], self.im_2_2)

        intervals = {repr(i) for i in self.im_2._query_interval((1, 7))}
        self.assertEqual(len(intervals), 2)
        self.assertSetEqual(intervals,
                            {repr(i) for i in self.im_2._intervals})

    def test_query_interval_upper_bound(self):
        intervals = list(self.im_2._query_interval((self.upper_bound-1,
                                                    self.upper_bound)))
        self.assertEqual(intervals, [self.im_2_1])

    def test_query(self):
        intervals = list(self.im_2.query(bounds=[(1, 5)],
                                         metadata={'gene': 'sagA'}))
        self.assertEqual(len(intervals), 1)
        self.assertEqual(intervals[0], self.im_2_1)

    def test_query_empty(self):
        intervals = list(self.im_1.query())
        self.assertEqual(len(intervals), 1)

    def test_query_no_hits(self):
        intervals = list(self.im_2.query(bounds=[(self.upper_bound, 200)]))
        self.assertEqual(len(intervals), 0)

        intervals = list(self.im_2.query(metadata={'gene': 'sagC'}))
        self.assertEqual(len(intervals), 0)

        intervals = list(self.im_2.query(bounds=[(1, 2)],
                                         metadata={'gene': 'sagC'}))
        self.assertEqual(len(intervals), 0)

    def test_query_interval_only(self):
        for loc in [[(1, 7)],
                    [(1, 2), (3, 4)]]:
            intervals = list(self.im_2.query(bounds=loc))
            self.assertEqual(len(intervals), 2)
            self.assertEqual(intervals[0], self.im_2_1)
            self.assertEqual(intervals[1], self.im_2_2)

    def test_query_metadata_only(self):
        intervals = list(self.im_2.query(metadata={'gene': 'sagB'}))
        self.assertEqual(len(intervals), 1)
        self.assertEqual(intervals[0], self.im_2_2)

        intervals = list(self.im_2.query(metadata={'bound': 0}))
        self.assertEqual(len(intervals), 2)
        self.assertEqual(intervals[0], self.im_2_1)
        self.assertEqual(intervals[1], self.im_2_2)

    def test_drop(self):
        intvl = self.im_2._intervals[0]
        self.im_2.drop([intvl])
        self.assertEqual(len(self.im_2._intervals), 1)
        self.assertEqual(self.im_2._intervals[0], self.im_2_2)
        # test the intvl was set to dropped
        self.assertTrue(intvl.dropped)

    def test_drop_all(self):
        self.im_2.drop(self.im_2._intervals)
        self.assertEqual(self.im_2, self.im_empty)

    def test_drop_negate(self):
        intvl = self.im_2._intervals[0]
        self.im_2.drop([intvl], negate=True)
        self.assertEqual(len(self.im_2._intervals), 1)
        self.assertEqual(self.im_2._intervals[0], intvl)
        # test the dropped intvl was set to dropped
        self.assertTrue(self.im_2_2.dropped)

    def test_reverse(self):
        self.im_2._reverse()
        Interval(
            interval_metadata=self.im_empty,
            bounds=[(0, 6), (8, 9)],
            metadata={'gene': 'sagA',  'bound': 0})
        Interval(
            interval_metadata=self.im_empty,
            bounds=[(5, 7)],
            metadata={'gene': 'sagB', 'bound': 0, 'spam': [0]})
        self.assertEqual(self.im_2, self.im_empty)

    def test_eq_ne(self):
        im1 = IntervalMetadata(10)
        im1.add(metadata={'gene': 'sagA', 'bound': '0'},
                bounds=[(0, 2), (4, 7)])
        im1.add(metadata={'gene': 'sagB', 'bound': '3'},
                bounds=[(3, 5)])

        # The ordering shouldn't matter
        im2 = IntervalMetadata(10)
        im2.add(metadata={'gene': 'sagB', 'bound': '3'},
                bounds=[(3, 5)])
        im2.add(metadata={'gene': 'sagA', 'bound': '0'},
                bounds=[(0, 2), (4, 7)])

        im3 = IntervalMetadata(10)
        im3.add(metadata={'gene': 'sagA', 'bound': '3'},
                bounds=[(0, 2), (4, 7)])
        im3.add(metadata={'gene': 'sagB', 'bound': '3'},
                bounds=[(3, 5)])

        self.assertReallyEqual(im1, im2)
        self.assertReallyNotEqual(im1, im3)

    def test_ne_diff_bounds(self):
        im1 = IntervalMetadata(10)
        im2 = IntervalMetadata(9)
        intvl = {'bounds': [(0, 1)], 'metadata': {'spam': 'foo'}}
        im1.add(**intvl)
        im2.add(**intvl)
        self.assertReallyNotEqual(im1, im2)

    def test_repr(self):
        exp = '''0 interval features
-------------------'''
        self.assertEqual(repr(self.im_empty), exp)

        self.im_empty.add([(1, 2)], metadata={'gene': 'sagA'})

        exp = ("1 interval feature\n"
               "------------------\n"
               r"Interval\(interval_metadata=<[0-9]+>, bounds=\[\(1, 2\)\], "
               r"fuzzy=\[\(False, False\)\], metadata={'gene': 'sagA'}\)")
        self.assertRegex(repr(self.im_empty), exp)

        self.im_empty.add([(3, 4)], metadata={'gene': 'sagB'})
        self.im_empty.add([(3, 4)], metadata={'gene': 'sagC'})
        self.im_empty.add([(3, 4)], metadata={'gene': 'sagD'})
        self.im_empty.add([(3, 4)], metadata={'gene': 'sagE'})
        self.im_empty.add([(3, 4)], metadata={'gene': 'sagF'})
        exp = ("6 interval features\n"
               "-------------------\n"
               r"Interval\(interval_metadata=<[0-9]+>, bounds=\[\(1, 2\)\], "
               r"fuzzy=\[\(False, False\)\], metadata={'gene': 'sagA'}\)\n"
               r"Interval\(interval_metadata=<[0-9]+>, bounds=\[\(3, 4\)\], "
               r"fuzzy=\[\(False, False\)\], metadata={'gene': 'sagB'}\)\n"
               r"...\n"
               r"Interval\(interval_metadata=<[0-9]+>, bounds=\[\(3, 4\)\], "
               r"fuzzy=\[\(False, False\)\], metadata={'gene': 'sagE'}\)\n"
               r"Interval\(interval_metadata=<[0-9]+>, bounds=\[\(3, 4\)\], "
               r"fuzzy=\[\(False, False\)\], metadata={'gene': 'sagF'}\)")
        self.assertRegex(repr(self.im_empty), exp)


if __name__ == '__main__':
    unittest.main()
