# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

import matplotlib as mpl
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import scipy.spatial.distance
from IPython.core.display import Image, SVG

import skbio.sequence.distance
from skbio import DistanceMatrix, Sequence
from skbio.stats.distance import (
    DissimilarityMatrixError, DistanceMatrixError, MissingIDError,
    DissimilarityMatrix, randdm)
from skbio.stats.distance._base import (_preprocess_input,
                                        _run_monte_carlo_stats)
from skbio.util import assert_data_frame_almost_equal
from skbio.util._testing import assert_series_almost_equal


class DissimilarityMatrixTestData(TestCase):
    def setUp(self):
        self.dm_1x1_data = [[0.0]]
        self.dm_2x2_data = [[0.0, 0.123], [0.123, 0.0]]
        self.dm_2x2_asym_data = [[0.0, 1.0], [-2.0, 0.0]]
        self.dm_3x3_data = [[0.0, 0.01, 4.2], [0.01, 0.0, 12.0],
                            [4.2, 12.0, 0.0]]
        self.dm_5x5_data = [[0, 1, 2, 3, 4],
                            [5, 0, 6, 7, 8],
                            [9, 1, 0, 2, 3],
                            [4, 5, 6, 0, 7],
                            [8, 9, 1, 2, 0]]


class DissimilarityMatrixTests(DissimilarityMatrixTestData):
    def setUp(self):
        super(DissimilarityMatrixTests, self).setUp()

        self.dm_1x1 = DissimilarityMatrix(self.dm_1x1_data, ['a'])
        self.dm_2x2 = DissimilarityMatrix(self.dm_2x2_data, ['a', 'b'])
        self.dm_2x2_asym = DissimilarityMatrix(self.dm_2x2_asym_data,
                                               ['a', 'b'])
        self.dm_3x3 = DissimilarityMatrix(self.dm_3x3_data, ['a', 'b', 'c'])
        self.dm_5x5 = DissimilarityMatrix(self.dm_5x5_data, list('abcde'))

        self.dms = [self.dm_1x1, self.dm_2x2, self.dm_2x2_asym, self.dm_3x3]
        self.dm_shapes = [(1, 1), (2, 2), (2, 2), (3, 3)]
        self.dm_sizes = [1, 4, 4, 9]
        self.dm_transposes = [
            self.dm_1x1, self.dm_2x2,
            DissimilarityMatrix([[0, -2], [1, 0]], ['a', 'b']), self.dm_3x3]
        self.dm_redundant_forms = [np.array(self.dm_1x1_data),
                                   np.array(self.dm_2x2_data),
                                   np.array(self.dm_2x2_asym_data),
                                   np.array(self.dm_3x3_data)]

    def test_avoid_copy_on_construction(self):
        # ((data, expect_copy))
        tests = (([[0, 1], [1, 0]], True),
                 ([(0, 1), (1, 0)], True),
                 (((0, 1), (1, 0)), True),
                 (np.array([[0, 1], [1, 0]], dtype='int'), True),
                 (np.array([[0, 1], [1, 0]], dtype='float'), False),
                 (np.array([[0, 1], [1, 0]], dtype=np.float32), False),
                 (np.array([[0, 1], [1, 0]], dtype=np.float64), False),
                 (np.array([[0, 1], [1, 0]], dtype='double'), False))

        for data, expect in tests:
            obj = DissimilarityMatrix(data)
            self.assertEqual(id(obj.data) != id(data), expect)

    def test_within(self):
        exp = pd.DataFrame([['a', 'a', 0.0],
                            ['a', 'c', 4.2],
                            ['c', 'a', 4.2],
                            ['c', 'c', 0.0]],
                           columns=['i', 'j', 'value'])
        obs = self.dm_3x3.within(['a', 'c'])
        pdt.assert_frame_equal(obs, exp)

    def test_within_order_stability(self):
        exp = pd.DataFrame([['a', 'a', 0.0],
                            ['a', 'c', 4.2],
                            ['c', 'a', 4.2],
                            ['c', 'c', 0.0]],
                           columns=['i', 'j', 'value'])

        # NOTE: order was changed from ['a', 'c'] to ['c', 'a']
        # but the output order in exp is consistent with
        # test_within
        obs = self.dm_3x3.within(['c', 'a'])
        pdt.assert_frame_equal(obs, exp)
        obs = self.dm_3x3.within(['a', 'c'])
        pdt.assert_frame_equal(obs, exp)

    def test_within_missing_id(self):
        with self.assertRaisesRegex(MissingIDError, "not found."):
            self.dm_3x3.within(['x', 'a'])

    def test_between(self):
        exp = pd.DataFrame([['b', 'a', 5.],
                            ['b', 'c', 6.],
                            ['b', 'e', 8.],
                            ['d', 'a', 4.],
                            ['d', 'c', 6.],
                            ['d', 'e', 7.]],
                           columns=['i', 'j', 'value'])

        obs = self.dm_5x5.between(['b', 'd'], ['a', 'c', 'e'])
        pdt.assert_frame_equal(obs, exp)

    def test_between_order_stability(self):
        exp = pd.DataFrame([['b', 'a', 5.],
                            ['b', 'c', 6.],
                            ['b', 'e', 8.],
                            ['d', 'a', 4.],
                            ['d', 'c', 6.],
                            ['d', 'e', 7.]],
                           columns=['i', 'j', 'value'])

        # varying the order of the "i" values, result remains consistent
        # with the test_between result
        obs = self.dm_5x5.between(['d', 'b'], ['a', 'c', 'e'])
        pdt.assert_frame_equal(obs, exp)

        # varying the order of the "j" values, result remains consistent
        # with the test_between result
        obs = self.dm_5x5.between(['b', 'd'], ['a', 'e', 'c'])
        pdt.assert_frame_equal(obs, exp)

        # varying the order of the "i" and "j" values, result remains
        # consistent with the test_between result
        obs = self.dm_5x5.between(['d', 'b'], ['a', 'e', 'c'])
        pdt.assert_frame_equal(obs, exp)

    def test_between_overlap(self):
        exp = pd.DataFrame([['b', 'a', 5.],
                            ['b', 'd', 7.],
                            ['b', 'e', 8.],
                            ['d', 'a', 4.],
                            ['d', 'd', 0.],
                            ['d', 'e', 7.]],
                           columns=['i', 'j', 'value'])

        # 'd' in i and j overlap
        with self.assertRaisesRegex(KeyError, ("This constraint can "
                                               "removed with "
                                               "allow_overlap=True.")):
            self.dm_5x5.between(['b', 'd'], ['a', 'd', 'e'])

        obs = self.dm_5x5.between(['b', 'd'], ['a', 'd', 'e'],
                                  allow_overlap=True)
        pdt.assert_frame_equal(obs, exp)

    def test_between_missing_id(self):
        with self.assertRaisesRegex(MissingIDError, "not found."):
            self.dm_3x3.between(['x', 'a'], ['a', 'b', 'c'])

        with self.assertRaisesRegex(MissingIDError, "not found."):
            self.dm_3x3.between(['a', 'b'], ['a', 'x', 'c'])

        with self.assertRaisesRegex(MissingIDError, "not found."):
            self.dm_3x3.between(['a', 'y'], ['a', 'x', 'c'])

    def test_stable_order(self):
        exp = np.array([1, 3, 4], dtype=int)
        obs = self.dm_5x5._stable_order(['d', 'e', 'b'])

        npt.assert_equal(obs, exp)

    def test_subset_to_dataframe(self):
        exp = pd.DataFrame([['b', 'a', 5.],
                            ['b', 'd', 7.],
                            ['b', 'e', 8.],
                            ['d', 'a', 4.],
                            ['d', 'd', 0.],
                            ['d', 'e', 7.]],
                           columns=['i', 'j', 'value'])

        obs = self.dm_5x5._subset_to_dataframe(['b', 'd'], ['a', 'd', 'e'])
        pdt.assert_frame_equal(obs, exp)

        # and the empty edge cases
        exp = pd.DataFrame([],
                           columns=['i', 'j', 'value'],
                           index=pd.RangeIndex(start=0, stop=0))

        obs = self.dm_5x5._subset_to_dataframe([], ['a', 'd', 'e'])
        pdt.assert_frame_equal(obs, exp, check_dtype=False)
        obs = self.dm_5x5._subset_to_dataframe(['b', 'd'], [])
        pdt.assert_frame_equal(obs, exp, check_dtype=False)
        obs = self.dm_5x5._subset_to_dataframe([], [])
        pdt.assert_frame_equal(obs, exp, check_dtype=False)

    def test_init_from_dm(self):
        ids = ['foo', 'bar', 'baz']

        # DissimilarityMatrix -> DissimilarityMatrix
        exp = DissimilarityMatrix(self.dm_3x3_data, ids)
        obs = DissimilarityMatrix(self.dm_3x3, ids)
        self.assertEqual(obs, exp)
        # Test that copy of data is not made.
        self.assertTrue(obs.data is self.dm_3x3.data)
        obs.data[0, 1] = 424242
        self.assertTrue(np.array_equal(obs.data, self.dm_3x3.data))

        # DistanceMatrix -> DissimilarityMatrix
        exp = DissimilarityMatrix(self.dm_3x3_data, ids)
        obs = DissimilarityMatrix(
            DistanceMatrix(self.dm_3x3_data, ('a', 'b', 'c')), ids)
        self.assertEqual(obs, exp)

        # DissimilarityMatrix -> DistanceMatrix
        with self.assertRaises(DistanceMatrixError):
            DistanceMatrix(self.dm_2x2_asym, ['foo', 'bar'])

    def test_init_non_hollow_dm(self):
        data = [[1, 1], [1, 1]]
        obs = DissimilarityMatrix(data, ['a', 'b'])
        self.assertTrue(np.array_equal(obs.data, data))

    def test_init_no_ids(self):
        exp = DissimilarityMatrix(self.dm_3x3_data, ('0', '1', '2'))
        obs = DissimilarityMatrix(self.dm_3x3_data)
        self.assertEqual(obs, exp)
        self.assertEqual(obs['1', '2'], 12.0)

    def test_init_invalid_input(self):
        # Empty data.
        with self.assertRaises(DissimilarityMatrixError):
            DissimilarityMatrix([], [])

        # Another type of empty data.
        with self.assertRaises(DissimilarityMatrixError):
            DissimilarityMatrix(np.empty((0, 0)), [])

        # Invalid number of dimensions.
        with self.assertRaises(DissimilarityMatrixError):
            DissimilarityMatrix([1, 2, 3], ['a'])

        # Dimensions don't match.
        with self.assertRaises(DissimilarityMatrixError):
            DissimilarityMatrix([[1, 2, 3]], ['a'])

        data = [[0, 1], [1, 0]]

        # Duplicate IDs.
        with self.assertRaises(DissimilarityMatrixError):
            DissimilarityMatrix(data, ['a', 'a'])

        # Number of IDs don't match dimensions.
        with self.assertRaises(DissimilarityMatrixError):
            DissimilarityMatrix(data, ['a', 'b', 'c'])
        with self.assertRaises(DissimilarityMatrixError):
            DissimilarityMatrix(data, [])

    def test_from_iterable_non_hollow_data(self):
        iterable = (x for x in range(4))

        exp = DissimilarityMatrix([[1, 1, 1, 1],
                                   [1, 1, 1, 1],
                                   [1, 1, 1, 1],
                                   [1, 1, 1, 1]])
        res = DissimilarityMatrix.from_iterable(iterable, lambda a, b: 1)
        self.assertEqual(res, exp)

    def test_from_iterable_asymmetric_data(self):
        iterable = (x for x in range(4))

        exp = DissimilarityMatrix([[0, 1, 2, 3],
                                   [-1, 0, 1, 2],
                                   [-2, -1, 0, 1],
                                   [-3, -2, -1, 0]])
        res = DissimilarityMatrix.from_iterable(iterable, lambda a, b: b - a)
        self.assertEqual(res, exp)

    def test_from_iterable_no_key(self):
        iterable = (x for x in range(4))

        exp = DissimilarityMatrix([[0, 1, 2, 3],
                                   [1, 0, 1, 2],
                                   [2, 1, 0, 1],
                                   [3, 2, 1, 0]])
        res = DissimilarityMatrix.from_iterable(iterable,
                                                lambda a, b: abs(b - a))
        self.assertEqual(res, exp)

    def test_from_iterable_with_key(self):
        iterable = (x for x in range(4))

        exp = DissimilarityMatrix([[0, 1, 2, 3],
                                   [1, 0, 1, 2],
                                   [2, 1, 0, 1],
                                   [3, 2, 1, 0]], ['0', '1', '4', '9'])
        res = DissimilarityMatrix.from_iterable(iterable,
                                                lambda a, b: abs(b - a),
                                                key=lambda x: str(x ** 2))
        self.assertEqual(res, exp)

    def test_from_iterable_empty(self):
        with self.assertRaises(DissimilarityMatrixError):
            DissimilarityMatrix.from_iterable([], lambda x: x)

    def test_from_iterable_single(self):
        exp = DissimilarityMatrix([[100]])
        res = DissimilarityMatrix.from_iterable(["boo"], lambda a, b: 100)
        self.assertEqual(res, exp)

    def test_from_iterable_with_keys(self):
        iterable = (x for x in range(4))

        exp = DissimilarityMatrix([[0, 1, 2, 3],
                                   [1, 0, 1, 2],
                                   [2, 1, 0, 1],
                                   [3, 2, 1, 0]], ['0', '1', '4', '9'])
        res = DissimilarityMatrix.from_iterable(iterable,
                                                lambda a, b: abs(b - a),
                                                keys=iter(['0', '1', '4', '9'])
                                                )
        self.assertEqual(res, exp)

    def test_from_iterable_with_key_and_keys(self):
        iterable = (x for x in range(4))
        with self.assertRaises(ValueError):
            DissimilarityMatrix.from_iterable(iterable,
                                              lambda a, b: abs(b - a),
                                              key=str,
                                              keys=['1', '2', '3', '4'])

    def test_from_iterable_scipy_hamming_metric_with_metadata(self):
        # test for #1254
        seqs = [
            Sequence('ACGT'),
            Sequence('ACGA', metadata={'id': 'seq1'}),
            Sequence('AAAA', metadata={'id': 'seq2'}),
            Sequence('AAAA', positional_metadata={'qual': range(4)})
        ]

        exp = DissimilarityMatrix([
            [0, 0.25, 0.75, 0.75],
            [0.25, 0.0, 0.5, 0.5],
            [0.75, 0.5, 0.0, 0.0],
            [0.75, 0.5, 0.0, 0.0]], ['a', 'b', 'c', 'd'])

        dm = DissimilarityMatrix.from_iterable(
            seqs,
            metric=scipy.spatial.distance.hamming,
            keys=['a', 'b', 'c', 'd'])

        self.assertEqual(dm, exp)

    def test_from_iterable_skbio_hamming_metric_with_metadata(self):
        # test for #1254
        seqs = [
            Sequence('ACGT'),
            Sequence('ACGA', metadata={'id': 'seq1'}),
            Sequence('AAAA', metadata={'id': 'seq2'}),
            Sequence('AAAA', positional_metadata={'qual': range(4)})
        ]

        exp = DissimilarityMatrix([
            [0, 0.25, 0.75, 0.75],
            [0.25, 0.0, 0.5, 0.5],
            [0.75, 0.5, 0.0, 0.0],
            [0.75, 0.5, 0.0, 0.0]], ['a', 'b', 'c', 'd'])

        dm = DissimilarityMatrix.from_iterable(
            seqs,
            metric=skbio.sequence.distance.hamming,
            keys=['a', 'b', 'c', 'd'])

        self.assertEqual(dm, exp)

    def test_data(self):
        for dm, exp in zip(self.dms, self.dm_redundant_forms):
            obs = dm.data
            self.assertTrue(np.array_equal(obs, exp))

        with self.assertRaises(AttributeError):
            self.dm_3x3.data = 'foo'

    def test_ids(self):
        obs = self.dm_3x3.ids
        self.assertEqual(obs, ('a', 'b', 'c'))

        # Test that we overwrite the existing IDs and that the ID index is
        # correctly rebuilt.
        new_ids = ['foo', 'bar', 'baz']
        self.dm_3x3.ids = new_ids
        obs = self.dm_3x3.ids
        self.assertEqual(obs, tuple(new_ids))
        self.assertTrue(np.array_equal(self.dm_3x3['bar'],
                                       np.array([0.01, 0.0, 12.0])))
        with self.assertRaises(MissingIDError):
            self.dm_3x3['b']

    def test_ids_invalid_input(self):
        with self.assertRaises(DissimilarityMatrixError):
            self.dm_3x3.ids = ['foo', 'bar']
        # Make sure that we can still use the dissimilarity matrix after trying
        # to be evil.
        obs = self.dm_3x3.ids
        self.assertEqual(obs, ('a', 'b', 'c'))

    def test_dtype(self):
        for dm in self.dms:
            self.assertEqual(dm.dtype, np.float64)

    def test_shape(self):
        for dm, shape in zip(self.dms, self.dm_shapes):
            self.assertEqual(dm.shape, shape)

    def test_size(self):
        for dm, size in zip(self.dms, self.dm_sizes):
            self.assertEqual(dm.size, size)

    def test_transpose(self):
        for dm, transpose in zip(self.dms, self.dm_transposes):
            self.assertEqual(dm.T, transpose)
            self.assertEqual(dm.transpose(), transpose)
            # We should get a reference to a different object back, even if the
            # transpose is the same as the original.
            self.assertTrue(dm.transpose() is not dm)

    def test_index(self):
        self.assertEqual(self.dm_3x3.index('a'), 0)
        self.assertEqual(self.dm_3x3.index('b'), 1)
        self.assertEqual(self.dm_3x3.index('c'), 2)

        with self.assertRaises(MissingIDError):
            self.dm_3x3.index('d')

        with self.assertRaises(MissingIDError):
            self.dm_3x3.index(1)

    def test_redundant_form(self):
        for dm, redundant in zip(self.dms, self.dm_redundant_forms):
            obs = dm.redundant_form()
            self.assertTrue(np.array_equal(obs, redundant))

    def test_copy(self):
        copy = self.dm_2x2.copy()
        self.assertEqual(copy, self.dm_2x2)
        self.assertFalse(copy.data is self.dm_2x2.data)
        # deepcopy doesn't actually create a copy of the IDs because it is a
        # tuple of strings, which is fully immutable.
        self.assertTrue(copy.ids is self.dm_2x2.ids)

        new_ids = ['hello', 'world']
        copy.ids = new_ids
        self.assertNotEqual(copy.ids, self.dm_2x2.ids)

        copy = self.dm_2x2.copy()
        copy.data[0, 1] = 0.0001
        self.assertFalse(np.array_equal(copy.data, self.dm_2x2.data))

    def test_filter_no_filtering(self):
        # Don't actually filter anything -- ensure we get back a different
        # object.
        obs = self.dm_3x3.filter(['a', 'b', 'c'])
        self.assertEqual(obs, self.dm_3x3)
        self.assertFalse(obs is self.dm_3x3)

    def test_filter_reorder(self):
        # Don't filter anything, but reorder the distance matrix.
        order = ['c', 'a', 'b']
        exp = DissimilarityMatrix(
            [[0, 4.2, 12], [4.2, 0, 0.01], [12, 0.01, 0]], order)
        obs = self.dm_3x3.filter(order)
        self.assertEqual(obs, exp)

    def test_filter_single_id(self):
        ids = ['b']
        exp = DissimilarityMatrix([[0]], ids)
        obs = self.dm_2x2_asym.filter(ids)
        self.assertEqual(obs, exp)

    def test_filter_asymmetric(self):
        # 2x2
        ids = ['b', 'a']
        exp = DissimilarityMatrix([[0, -2], [1, 0]], ids)
        obs = self.dm_2x2_asym.filter(ids)
        self.assertEqual(obs, exp)

        # 3x3
        dm = DissimilarityMatrix([[0, 10, 53], [42, 0, 22.5], [53, 1, 0]],
                                 ('bro', 'brah', 'breh'))
        ids = ['breh', 'brah']
        exp = DissimilarityMatrix([[0, 1], [22.5, 0]], ids)
        obs = dm.filter(ids)
        self.assertEqual(obs, exp)

    def test_filter_subset(self):
        ids = ('c', 'a')
        exp = DissimilarityMatrix([[0, 4.2], [4.2, 0]], ids)
        obs = self.dm_3x3.filter(ids)
        self.assertEqual(obs, exp)

        ids = ('b', 'a')
        exp = DissimilarityMatrix([[0, 0.01], [0.01, 0]], ids)
        obs = self.dm_3x3.filter(ids)
        self.assertEqual(obs, exp)

        # 4x4
        dm = DissimilarityMatrix([[0, 1, 55, 7], [1, 0, 16, 1],
                                  [55, 16, 0, 23], [7, 1, 23, 0]])
        ids = np.asarray(['3', '0', '1'])
        exp = DissimilarityMatrix([[0, 7, 1], [7, 0, 1], [1, 1, 0]], ids)
        obs = dm.filter(ids)
        self.assertEqual(obs, exp)

    def test_filter_duplicate_ids(self):
        with self.assertRaises(DissimilarityMatrixError):
            self.dm_3x3.filter(['c', 'a', 'c'])

    def test_filter_missing_ids(self):
        with self.assertRaises(MissingIDError):
            self.dm_3x3.filter(['c', 'bro'])

    def test_filter_missing_ids_strict_false(self):
        # no exception should be raised
        ids = ('c', 'a')
        exp = DissimilarityMatrix([[0, 4.2], [4.2, 0]], ids)
        obs = self.dm_3x3.filter(['c', 'a', 'not found'], strict=False)
        self.assertEqual(obs, exp)

    def test_filter_empty_ids(self):
        with self.assertRaises(DissimilarityMatrixError):
            self.dm_3x3.filter([])

    def test_plot_default(self):
        fig = self.dm_1x1.plot()
        self.assertIsInstance(fig, mpl.figure.Figure)
        axes = fig.get_axes()
        self.assertEqual(len(axes), 2)
        ax = axes[0]
        self.assertEqual(ax.get_title(), '')
        xticks = []
        for tick in ax.get_xticklabels():
            xticks.append(tick.get_text())
        self.assertEqual(xticks, ['a'])
        yticks = []
        for tick in ax.get_yticklabels():
            yticks.append(tick.get_text())
        self.assertEqual(yticks, ['a'])

    def test_plot_no_default(self):
        ids = ['0', 'one', '2', 'three', '4.000']
        data = ([0, 1, 2, 3, 4], [1, 0, 1, 2, 3], [2, 1, 0, 1, 2],
                [3, 2, 1, 0, 1], [4, 3, 2, 1, 0])
        dm = DissimilarityMatrix(data, ids)
        fig = dm.plot(cmap='Reds', title='Testplot')
        self.assertIsInstance(fig, mpl.figure.Figure)
        axes = fig.get_axes()
        self.assertEqual(len(axes), 2)
        ax = axes[0]
        self.assertEqual(ax.get_title(), 'Testplot')
        xticks = []
        for tick in ax.get_xticklabels():
            xticks.append(tick.get_text())
        self.assertEqual(xticks, ['0', 'one', '2', 'three', '4.000'])
        yticks = []
        for tick in ax.get_yticklabels():
            yticks.append(tick.get_text())
        self.assertEqual(yticks, ['0', 'one', '2', 'three', '4.000'])

    def test_repr_png(self):
        dm = self.dm_1x1
        obs = dm._repr_png_()
        self.assertIsInstance(obs, bytes)
        self.assertTrue(len(obs) > 0)

    def test_repr_svg(self):
        obs = self.dm_1x1._repr_svg_()
        self.assertIsInstance(obs, str)
        self.assertTrue(len(obs) > 0)

    def test_png(self):
        dm = self.dm_1x1
        self.assertIsInstance(dm.png, Image)

    def test_svg(self):
        dm = self.dm_1x1
        self.assertIsInstance(dm.svg, SVG)

    def test_to_data_frame_1x1(self):
        df = self.dm_1x1.to_data_frame()
        exp = pd.DataFrame([[0.0]], index=['a'], columns=['a'])
        assert_data_frame_almost_equal(df, exp)

    def test_to_data_frame_3x3(self):
        df = self.dm_3x3.to_data_frame()
        exp = pd.DataFrame([[0.0, 0.01, 4.2],
                            [0.01, 0.0, 12.0],
                            [4.2, 12.0, 0.0]],
                           index=['a', 'b', 'c'], columns=['a', 'b', 'c'])
        assert_data_frame_almost_equal(df, exp)

    def test_to_data_frame_default_ids(self):
        df = DissimilarityMatrix(self.dm_2x2_data).to_data_frame()
        exp = pd.DataFrame([[0.0, 0.123],
                            [0.123, 0.0]],
                           index=['0', '1'], columns=['0', '1'])
        assert_data_frame_almost_equal(df, exp)

    def test_str(self):
        for dm in self.dms:
            obs = str(dm)
            # Do some very light testing here to make sure we're getting a
            # non-empty string back. We don't want to test the exact
            # formatting.
            self.assertTrue(obs)

    def test_eq(self):
        for dm in self.dms:
            copy = dm.copy()
            self.assertTrue(dm == dm)
            self.assertTrue(copy == copy)
            self.assertTrue(dm == copy)
            self.assertTrue(copy == dm)

        self.assertFalse(self.dm_1x1 == self.dm_3x3)

    def test_ne(self):
        # Wrong class.
        self.assertTrue(self.dm_3x3 != 'foo')

        # Wrong shape.
        self.assertTrue(self.dm_3x3 != self.dm_1x1)

        # Wrong IDs.
        other = self.dm_3x3.copy()
        other.ids = ['foo', 'bar', 'baz']
        self.assertTrue(self.dm_3x3 != other)

        # Wrong data.
        other = self.dm_3x3.copy()
        other.data[1, 0] = 42.42
        self.assertTrue(self.dm_3x3 != other)

        self.assertFalse(self.dm_2x2 != self.dm_2x2)

    def test_contains(self):
        self.assertTrue('a' in self.dm_3x3)
        self.assertTrue('b' in self.dm_3x3)
        self.assertTrue('c' in self.dm_3x3)
        self.assertFalse('d' in self.dm_3x3)

    def test_getslice(self):
        # Slice of first dimension only. Test that __getslice__ defers to
        # __getitem__.
        obs = self.dm_2x2[1:]
        self.assertTrue(np.array_equal(obs, np.array([[0.123, 0.0]])))
        self.assertEqual(type(obs), np.ndarray)

    def test_getitem_by_id(self):
        obs = self.dm_1x1['a']
        self.assertTrue(np.array_equal(obs, np.array([0.0])))

        obs = self.dm_2x2_asym['b']
        self.assertTrue(np.array_equal(obs, np.array([-2.0, 0.0])))

        obs = self.dm_3x3['c']
        self.assertTrue(np.array_equal(obs, np.array([4.2, 12.0, 0.0])))

        with self.assertRaises(MissingIDError):
            self.dm_2x2['c']

    def test_getitem_by_id_pair(self):
        # Same object.
        self.assertEqual(self.dm_1x1['a', 'a'], 0.0)

        # Different objects (symmetric).
        self.assertEqual(self.dm_3x3['b', 'c'], 12.0)
        self.assertEqual(self.dm_3x3['c', 'b'], 12.0)

        # Different objects (asymmetric).
        self.assertEqual(self.dm_2x2_asym['a', 'b'], 1.0)
        self.assertEqual(self.dm_2x2_asym['b', 'a'], -2.0)

        with self.assertRaises(MissingIDError):
            self.dm_2x2['a', 'c']

    def test_getitem_ndarray_indexing(self):
        # Single element access.
        obs = self.dm_3x3[0, 1]
        self.assertEqual(obs, 0.01)

        # Single element access (via two __getitem__ calls).
        obs = self.dm_3x3[0][1]
        self.assertEqual(obs, 0.01)

        # Row access.
        obs = self.dm_3x3[1]
        self.assertTrue(np.array_equal(obs, np.array([0.01, 0.0, 12.0])))
        self.assertEqual(type(obs), np.ndarray)

        # Grab all data.
        obs = self.dm_3x3[:, :]
        self.assertTrue(np.array_equal(obs, self.dm_3x3.data))
        self.assertEqual(type(obs), np.ndarray)

        with self.assertRaises(IndexError):
            self.dm_3x3[:, 3]

    def test_validate_invalid_dtype(self):
        with self.assertRaises(DissimilarityMatrixError):
            self.dm_3x3._validate(np.array([[0, 42], [42, 0]]), ['a', 'b'])


class DistanceMatrixTests(DissimilarityMatrixTestData):
    def setUp(self):
        super(DistanceMatrixTests, self).setUp()

        self.dm_1x1 = DistanceMatrix(self.dm_1x1_data, ['a'])
        self.dm_2x2 = DistanceMatrix(self.dm_2x2_data, ['a', 'b'])
        self.dm_3x3 = DistanceMatrix(self.dm_3x3_data, ['a', 'b', 'c'])

        self.dms = [self.dm_1x1, self.dm_2x2, self.dm_3x3]
        self.dm_condensed_forms = [np.array([]), np.array([0.123]),
                                   np.array([0.01, 4.2, 12.0])]

    def test_init_from_condensed_form(self):
        data = [1, 2, 3]
        exp = DistanceMatrix([[0, 1, 2],
                              [1, 0, 3],
                              [2, 3, 0]], ['0', '1', '2'])
        res = DistanceMatrix(data)
        self.assertEqual(exp, res)

    def test_init_invalid_input(self):
        # Asymmetric.
        data = [[0.0, 2.0], [1.0, 0.0]]
        with self.assertRaises(DistanceMatrixError):
            DistanceMatrix(data, ['a', 'b'])

        # Non-hollow
        data = [[1.0, 2.0], [2.0, 1.0]]
        with self.assertRaises(DistanceMatrixError):
            DistanceMatrix(data, ['a', 'b'])

        # Ensure that the superclass validation is still being performed.
        with self.assertRaises(DissimilarityMatrixError):
            DistanceMatrix([[1, 2, 3]], ['a'])

    def test_init_nans(self):
        with self.assertRaisesRegex(DistanceMatrixError, r'NaNs'):
            DistanceMatrix([[0.0, np.nan], [np.nan, 0.0]], ['a', 'b'])

    def test_from_iterable_no_key(self):
        iterable = (x for x in range(4))

        exp = DistanceMatrix([[0, 1, 2, 3],
                              [1, 0, 1, 2],
                              [2, 1, 0, 1],
                              [3, 2, 1, 0]])
        res = DistanceMatrix.from_iterable(iterable, lambda a, b: abs(b - a))
        self.assertEqual(res, exp)

    def test_from_iterable_validate_equal_valid_data(self):
        validate_true = DistanceMatrix.from_iterable((x for x in range(4)),
                                                     lambda a, b: abs(b - a),
                                                     validate=True)
        validate_false = DistanceMatrix.from_iterable((x for x in range(4)),
                                                      lambda a, b: abs(b - a),
                                                      validate=False)
        self.assertEqual(validate_true, validate_false)

    def test_from_iterable_validate_false(self):
        iterable = (x for x in range(4))

        exp = DistanceMatrix([[0, 1, 2, 3],
                              [1, 0, 1, 2],
                              [2, 1, 0, 1],
                              [3, 2, 1, 0]])
        res = DistanceMatrix.from_iterable(iterable, lambda a, b: abs(b - a),
                                           validate=False)
        self.assertEqual(res, exp)

    def test_from_iterable_validate_non_hollow(self):
        iterable = (x for x in range(4))
        with self.assertRaises(DistanceMatrixError):
            DistanceMatrix.from_iterable(iterable, lambda a, b: 1)

    def test_from_iterable_validate_false_non_symmetric(self):
        exp = DistanceMatrix([[0, 1, 2, 3],
                              [1, 0, 1, 2],
                              [2, 1, 0, 1],
                              [3, 2, 1, 0]])
        res = DistanceMatrix.from_iterable((x for x in range(4)),
                                           lambda a, b: a - b,
                                           validate=False)
        self.assertEqual(res, exp)

    def test_from_iterable_validate_asym(self):
        iterable = (x for x in range(4))
        with self.assertRaises(DistanceMatrixError):
            DistanceMatrix.from_iterable(iterable, lambda a, b: b - a)

    def test_from_iterable_with_key(self):
        iterable = (x for x in range(4))

        exp = DistanceMatrix([[0, 1, 2, 3],
                              [1, 0, 1, 2],
                              [2, 1, 0, 1],
                              [3, 2, 1, 0]], ['0', '1', '4', '9'])
        res = DistanceMatrix.from_iterable(iterable, lambda a, b: abs(b - a),
                                           key=lambda x: str(x**2))
        self.assertEqual(res, exp)

    def test_from_iterable_empty(self):
        with self.assertRaises(DissimilarityMatrixError):
            DistanceMatrix.from_iterable([], lambda x: x)

    def test_from_iterable_single(self):
        exp = DistanceMatrix([[0]])
        res = DistanceMatrix.from_iterable(["boo"], lambda a, b: 0)
        self.assertEqual(res, exp)

    def test_from_iterable_with_keys(self):
        iterable = (x for x in range(4))

        exp = DistanceMatrix([[0, 1, 2, 3],
                              [1, 0, 1, 2],
                              [2, 1, 0, 1],
                              [3, 2, 1, 0]], ['0', '1', '4', '9'])
        res = DistanceMatrix.from_iterable(iterable, lambda a, b: abs(b - a),
                                           keys=iter(['0', '1', '4', '9']))
        self.assertEqual(res, exp)

    def test_from_iterable_with_key_and_keys(self):
        iterable = (x for x in range(4))
        with self.assertRaises(ValueError):
            DistanceMatrix.from_iterable(iterable, lambda a, b: abs(b - a),
                                         key=str, keys=['1', '2', '3', '4'])

    def test_from_iterable_scipy_hamming_metric_with_metadata(self):
        # test for #1254
        seqs = [
            Sequence('ACGT'),
            Sequence('ACGA', metadata={'id': 'seq1'}),
            Sequence('AAAA', metadata={'id': 'seq2'}),
            Sequence('AAAA', positional_metadata={'qual': range(4)})
        ]

        exp = DistanceMatrix([
            [0, 0.25, 0.75, 0.75],
            [0.25, 0.0, 0.5, 0.5],
            [0.75, 0.5, 0.0, 0.0],
            [0.75, 0.5, 0.0, 0.0]], ['a', 'b', 'c', 'd'])

        dm = DistanceMatrix.from_iterable(
            seqs,
            metric=scipy.spatial.distance.hamming,
            keys=['a', 'b', 'c', 'd'])

        self.assertEqual(dm, exp)

    def test_from_iterable_skbio_hamming_metric_with_metadata(self):
        # test for #1254
        seqs = [
            Sequence('ACGT'),
            Sequence('ACGA', metadata={'id': 'seq1'}),
            Sequence('AAAA', metadata={'id': 'seq2'}),
            Sequence('AAAA', positional_metadata={'qual': range(4)})
        ]

        exp = DistanceMatrix([
            [0, 0.25, 0.75, 0.75],
            [0.25, 0.0, 0.5, 0.5],
            [0.75, 0.5, 0.0, 0.0],
            [0.75, 0.5, 0.0, 0.0]], ['a', 'b', 'c', 'd'])

        dm = DistanceMatrix.from_iterable(
            seqs,
            metric=skbio.sequence.distance.hamming,
            keys=['a', 'b', 'c', 'd'])

        self.assertEqual(dm, exp)

    def test_condensed_form(self):
        for dm, condensed in zip(self.dms, self.dm_condensed_forms):
            obs = dm.condensed_form()
            self.assertTrue(np.array_equal(obs, condensed))

    def test_permute_condensed(self):
        # Can't really permute a 1x1 or 2x2...
        for _ in range(2):
            obs = self.dm_1x1.permute(condensed=True)
            npt.assert_equal(obs, np.array([]))

        for _ in range(2):
            obs = self.dm_2x2.permute(condensed=True)
            npt.assert_equal(obs, np.array([0.123]))

        dm_copy = self.dm_3x3.copy()

        np.random.seed(0)

        obs = self.dm_3x3.permute(condensed=True)
        npt.assert_equal(obs, np.array([12.0, 4.2, 0.01]))

        obs = self.dm_3x3.permute(condensed=True)
        npt.assert_equal(obs, np.array([4.2, 12.0, 0.01]))

        # Ensure dm hasn't changed after calling permute() on it a couple of
        # times.
        self.assertEqual(self.dm_3x3, dm_copy)

    def test_permute_not_condensed(self):
        obs = self.dm_1x1.permute()
        self.assertEqual(obs, self.dm_1x1)
        self.assertFalse(obs is self.dm_1x1)

        obs = self.dm_2x2.permute()
        self.assertEqual(obs, self.dm_2x2)
        self.assertFalse(obs is self.dm_2x2)

        np.random.seed(0)

        exp = DistanceMatrix([[0, 12, 4.2],
                              [12, 0, 0.01],
                              [4.2, 0.01, 0]], self.dm_3x3.ids)
        obs = self.dm_3x3.permute()
        self.assertEqual(obs, exp)

        exp = DistanceMatrix([[0, 4.2, 12],
                              [4.2, 0, 0.01],
                              [12, 0.01, 0]], self.dm_3x3.ids)
        obs = self.dm_3x3.permute()
        self.assertEqual(obs, exp)

    def test_eq(self):
        # Compare DistanceMatrix to DissimilarityMatrix, where both have the
        # same data and IDs.
        eq_dm = DissimilarityMatrix(self.dm_3x3_data, ['a', 'b', 'c'])
        self.assertTrue(self.dm_3x3 == eq_dm)
        self.assertTrue(eq_dm == self.dm_3x3)

    def test_to_series_1x1(self):
        series = self.dm_1x1.to_series()

        exp = pd.Series([], index=[])
        assert_series_almost_equal(series, exp)

    def test_to_series_2x2(self):
        series = self.dm_2x2.to_series()

        exp = pd.Series([0.123], index=pd.Index([('a', 'b')]))
        assert_series_almost_equal(series, exp)

    def test_to_series_4x4(self):
        dm = DistanceMatrix([
            [0.0, 0.2, 0.3, 0.4],
            [0.2, 0.0, 0.5, 0.6],
            [0.3, 0.5, 0.0, 0.7],
            [0.4, 0.6, 0.7, 0.0]], ['a', 'b', 'c', 'd'])

        series = dm.to_series()

        exp = pd.Series([0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
                        index=pd.Index([('a', 'b'), ('a', 'c'), ('a', 'd'),
                                        ('b', 'c'), ('b', 'd'), ('c', 'd')]))
        assert_series_almost_equal(series, exp)

    def test_to_series_default_ids(self):
        series = DistanceMatrix(self.dm_2x2_data).to_series()

        exp = pd.Series([0.123], index=pd.Index([('0', '1')]))
        assert_series_almost_equal(series, exp)


class RandomDistanceMatrixTests(TestCase):
    def test_default_usage(self):
        exp = DistanceMatrix(np.asarray([[0.0]]), ['1'])
        obs = randdm(1)
        self.assertEqual(obs, exp)

        obs = randdm(2)
        self.assertEqual(obs.shape, (2, 2))
        self.assertEqual(obs.ids, ('1', '2'))

        obs1 = randdm(5)
        num_trials = 10
        found_diff = False
        for _ in range(num_trials):
            obs2 = randdm(5)

            if obs1 != obs2:
                found_diff = True
                break

        self.assertTrue(found_diff)

    def test_large_matrix_for_symmetry(self):
        obs3 = randdm(100)
        self.assertEqual(obs3, obs3.T)

    def test_ids(self):
        ids = ['foo', 'bar', 'baz']
        obs = randdm(3, ids=ids)
        self.assertEqual(obs.shape, (3, 3))
        self.assertEqual(obs.ids, tuple(ids))

    def test_constructor(self):
        exp = DissimilarityMatrix(np.asarray([[0.0]]), ['1'])
        obs = randdm(1, constructor=DissimilarityMatrix)
        self.assertEqual(obs, exp)
        self.assertEqual(type(obs), DissimilarityMatrix)

    def test_random_fn(self):
        def myrand(num_rows, num_cols):
            # One dm to rule them all...
            data = np.empty((num_rows, num_cols))
            data.fill(42)
            return data

        exp = DistanceMatrix(np.asarray([[0, 42, 42], [42, 0, 42],
                                         [42, 42, 0]]), ['1', '2', '3'])
        obs = randdm(3, random_fn=myrand)
        self.assertEqual(obs, exp)

    def test_invalid_input(self):
        # Invalid dimensions.
        with self.assertRaises(DissimilarityMatrixError):
            randdm(0)

        # Invalid dimensions.
        with self.assertRaises(ValueError):
            randdm(-1)

        # Invalid number of IDs.
        with self.assertRaises(DissimilarityMatrixError):
            randdm(2, ids=['foo'])


class CategoricalStatsHelperFunctionTests(TestCase):
    def setUp(self):
        self.dm = DistanceMatrix([[0.0, 1.0, 2.0],
                                  [1.0, 0.0, 3.0],
                                  [2.0, 3.0, 0.0]], ['a', 'b', 'c'])
        self.grouping = [1, 2, 1]
        # Ordering of IDs shouldn't matter, nor should extra IDs.
        self.df = pd.read_csv(
            io.StringIO('ID,Group\nb,Group2\na,Group1\nc,Group1\nd,Group3'),
            index_col=0)
        self.df_missing_id = pd.read_csv(
            io.StringIO('ID,Group\nb,Group2\nc,Group1'), index_col=0)

    def test_preprocess_input_with_valid_input(self):
        # Should obtain same result using grouping vector or data frame.
        exp = (3, 2, np.array([0, 1, 0]),
               (np.array([0, 0, 1]), np.array([1, 2, 2])),
               np.array([1., 2., 3.]))

        obs = _preprocess_input(self.dm, self.grouping, None)
        npt.assert_equal(obs, exp)

        obs = _preprocess_input(self.dm, self.df, 'Group')
        npt.assert_equal(obs, exp)

    def test_preprocess_input_raises_error(self):
        # Requires a DistanceMatrix.
        with self.assertRaises(TypeError):
            _preprocess_input(
                DissimilarityMatrix([[0, 2], [3, 0]], ['a', 'b']),
                [1, 2], None)

        # Requires column if DataFrame.
        with self.assertRaises(ValueError):
            _preprocess_input(self.dm, self.df, None)

        # Cannot provide column if not data frame.
        with self.assertRaises(ValueError):
            _preprocess_input(self.dm, self.grouping, 'Group')

        # Column must exist in data frame.
        with self.assertRaises(ValueError):
            _preprocess_input(self.dm, self.df, 'foo')

        # All distance matrix IDs must be in data frame.
        with self.assertRaises(ValueError):
            _preprocess_input(self.dm, self.df_missing_id, 'Group')

        # Grouping vector length must match number of objects in dm.
        with self.assertRaises(ValueError):
            _preprocess_input(self.dm, [1, 2], None)

        # Grouping vector cannot have only unique values.
        with self.assertRaises(ValueError):
            _preprocess_input(self.dm, [1, 2, 3], None)

        # Grouping vector cannot have only a single group.
        with self.assertRaises(ValueError):
            _preprocess_input(self.dm, [1, 1, 1], None)

    def test_run_monte_carlo_stats_with_permutations(self):
        obs = _run_monte_carlo_stats(lambda e: 42, self.grouping, 50)
        npt.assert_equal(obs, (42, 1.0))

    def test_run_monte_carlo_stats_no_permutations(self):
        obs = _run_monte_carlo_stats(lambda e: 42, self.grouping, 0)
        npt.assert_equal(obs, (42, np.nan))

    def test_run_monte_carlo_stats_invalid_permutations(self):
        with self.assertRaises(ValueError):
            _run_monte_carlo_stats(lambda e: 42, self.grouping, -1)


if __name__ == '__main__':
    main()
