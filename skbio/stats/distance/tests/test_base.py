# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main, skipUnless

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import scipy.spatial.distance
from scipy.sparse import csr_array, issparse

try:
    import matplotlib as mpl
except ImportError:
    has_matplotlib = False
else:
    has_matplotlib = True

import skbio.sequence.distance
from skbio import DistanceMatrix, Sequence
from skbio.stats.distance import (
    PairwiseMatrixError,
    DistanceMatrixError,
    SymmetricMatrixError,
    MissingIDError,
    PairwiseMatrix,
    SymmetricMatrix,
    randdm,
)
from skbio.stats.distance._base import _preprocess_input, _run_monte_carlo_stats
from skbio.stats.distance._utils import is_symmetric_and_hollow
from skbio.util import assert_data_frame_almost_equal
from skbio.util._testing import assert_series_almost_equal


class PairwiseMatrixTestData:
    def setUp(self):
        self.dm_1x1_data = [[0.0]]
        self.dm_2x2_data = [[0.0, 0.123], [0.123, 0.0]]
        self.dm_2x2_asym_data = [[0.0, 1.0], [-2.0, 0.0]]
        self.dm_3x3_data = [[0.0, 0.01, 4.2], [0.01, 0.0, 12.0], [4.2, 12.0, 0.0]]
        self.dm_5x5_data = [
            [0, 1, 2, 3, 4],
            [5, 0, 6, 7, 8],
            [9, 1, 0, 2, 3],
            [4, 5, 6, 0, 7],
            [8, 9, 1, 2, 0],
        ]


class PairwiseMatrixTestBase(PairwiseMatrixTestData):
    @classmethod
    def get_matrix_class(cls):
        return None

    def setUp(self):
        super(PairwiseMatrixTestBase, self).setUp()
        self.matobj = self.get_matrix_class()
        self.dm_1x1 = self.matobj(self.dm_1x1_data, ["a"])
        self.dm_2x2 = self.matobj(self.dm_2x2_data, ["a", "b"])
        self.dm_2x2_asym = self.matobj(self.dm_2x2_asym_data, ["a", "b"])
        self.dm_3x3 = self.matobj(self.dm_3x3_data, ["a", "b", "c"])
        self.dm_5x5 = self.matobj(self.dm_5x5_data, list("abcde"))

        self.dms = [self.dm_1x1, self.dm_2x2, self.dm_2x2_asym, self.dm_3x3]
        self.dm_shapes = [(1, 1), (2, 2), (2, 2), (3, 3)]
        self.dm_sizes = [1, 4, 4, 9]
        self.dm_transposes = [
            self.dm_1x1,
            self.dm_2x2,
            self.matobj([[0, -2], [1, 0]], ["a", "b"]),
            self.dm_3x3,
        ]
        self.dm_redundant_forms = [
            np.array(self.dm_1x1_data),
            np.array(self.dm_2x2_data),
            np.array(self.dm_2x2_asym_data),
            np.array(self.dm_3x3_data),
        ]

    def test_avoid_copy_on_construction(self):
        # ((data, expect_copy))
        tests = (
            ([[0, 1], [1, 0]], True),
            ([(0, 1), (1, 0)], True),
            (((0, 1), (1, 0)), True),
            (np.array([[0, 1], [1, 0]], dtype="int"), True),
            (np.array([[0, 1], [1, 0]], dtype="float"), False),
            (np.array([[0, 1], [1, 0]], dtype=np.float32), False),
            (np.array([[0, 1], [1, 0]], dtype=np.float64), False),
            (np.array([[0, 1], [1, 0]], dtype="double"), False),
        )

        for data, expect in tests:
            obj = PairwiseMatrix(data)
            self.assertEqual(id(obj.data) != id(data), expect)

    def test_within(self):
        exp = pd.DataFrame(
            [["a", "a", 0.0], ["a", "c", 4.2], ["c", "a", 4.2], ["c", "c", 0.0]],
            columns=["i", "j", "value"],
        )
        obs = self.dm_3x3.within(["a", "c"])
        pdt.assert_frame_equal(obs, exp)

    def test_within_order_stability(self):
        exp = pd.DataFrame(
            [["a", "a", 0.0], ["a", "c", 4.2], ["c", "a", 4.2], ["c", "c", 0.0]],
            columns=["i", "j", "value"],
        )

        # NOTE: order was changed from ['a', 'c'] to ['c', 'a']
        # but the output order in exp is consistent with
        # test_within
        obs = self.dm_3x3.within(["c", "a"])
        pdt.assert_frame_equal(obs, exp)
        obs = self.dm_3x3.within(["a", "c"])
        pdt.assert_frame_equal(obs, exp)

    def test_within_missing_id(self):
        with self.assertRaisesRegex(MissingIDError, "not found."):
            self.dm_3x3.within(["x", "a"])

    def test_between(self):
        exp = pd.DataFrame(
            [
                ["b", "a", 5.0],
                ["b", "c", 6.0],
                ["b", "e", 8.0],
                ["d", "a", 4.0],
                ["d", "c", 6.0],
                ["d", "e", 7.0],
            ],
            columns=["i", "j", "value"],
        )

        obs = self.dm_5x5.between(["b", "d"], ["a", "c", "e"])
        pdt.assert_frame_equal(obs, exp)

    def test_between_order_stability(self):
        exp = pd.DataFrame(
            [
                ["b", "a", 5.0],
                ["b", "c", 6.0],
                ["b", "e", 8.0],
                ["d", "a", 4.0],
                ["d", "c", 6.0],
                ["d", "e", 7.0],
            ],
            columns=["i", "j", "value"],
        )

        # varying the order of the "i" values, result remains consistent
        # with the test_between result
        obs = self.dm_5x5.between(["d", "b"], ["a", "c", "e"])
        pdt.assert_frame_equal(obs, exp)

        # varying the order of the "j" values, result remains consistent
        # with the test_between result
        obs = self.dm_5x5.between(["b", "d"], ["a", "e", "c"])
        pdt.assert_frame_equal(obs, exp)

        # varying the order of the "i" and "j" values, result remains
        # consistent with the test_between result
        obs = self.dm_5x5.between(["d", "b"], ["a", "e", "c"])
        pdt.assert_frame_equal(obs, exp)

    def test_between_overlap(self):
        exp = pd.DataFrame(
            [
                ["b", "a", 5.0],
                ["b", "d", 7.0],
                ["b", "e", 8.0],
                ["d", "a", 4.0],
                ["d", "d", 0.0],
                ["d", "e", 7.0],
            ],
            columns=["i", "j", "value"],
        )

        # 'd' in i and j overlap
        with self.assertRaisesRegex(
            KeyError, ("This constraint can removed with allow_overlap=True.")
        ):
            self.dm_5x5.between(["b", "d"], ["a", "d", "e"])

        obs = self.dm_5x5.between(["b", "d"], ["a", "d", "e"], allow_overlap=True)
        pdt.assert_frame_equal(obs, exp)

    def test_between_missing_id(self):
        with self.assertRaisesRegex(MissingIDError, "not found."):
            self.dm_3x3.between(["x", "a"], ["a", "b", "c"])

        with self.assertRaisesRegex(MissingIDError, "not found."):
            self.dm_3x3.between(["a", "b"], ["a", "x", "c"])

        with self.assertRaisesRegex(MissingIDError, "not found."):
            self.dm_3x3.between(["a", "y"], ["a", "x", "c"])

    def test_stable_order(self):
        exp = np.array([1, 3, 4], dtype=int)
        obs = self.dm_5x5._stable_order(["d", "e", "b"])

        npt.assert_equal(obs, exp)

    def test_matrix_from_matrix(self):
        # pairwise from pairwise
        dm = self.matobj(self.dm_5x5)
        npt.assert_equal(dm.data, self.dm_5x5.data)
        self.assertEqual(dm.ids, self.dm_5x5.ids)

        # pairwise from symmetric
        sm = SymmetricMatrix(self.dm_3x3)
        dm = self.matobj(sm)
        npt.assert_equal(dm.data, self.dm_3x3.data)
        self.assertEqual(dm.ids, self.dm_3x3.ids)

        # pairwise from distance
        dm = DistanceMatrix(self.dm_3x3)
        dm = self.matobj(dm)
        npt.assert_equal(dm.data, self.dm_3x3.data)
        self.assertEqual(dm.ids, self.dm_3x3.ids)

    def test_subset_to_dataframe(self):
        exp = pd.DataFrame(
            [
                ["b", "a", 5.0],
                ["b", "d", 7.0],
                ["b", "e", 8.0],
                ["d", "a", 4.0],
                ["d", "d", 0.0],
                ["d", "e", 7.0],
            ],
            columns=["i", "j", "value"],
        )

        obs = self.dm_5x5._subset_to_dataframe(["b", "d"], ["a", "d", "e"])
        pdt.assert_frame_equal(obs, exp)

        # and the empty edge cases
        exp = pd.DataFrame(
            [], columns=["i", "j", "value"], index=pd.RangeIndex(start=0, stop=0)
        )

        obs = self.dm_5x5._subset_to_dataframe([], ["a", "d", "e"])
        pdt.assert_frame_equal(obs, exp, check_dtype=False)
        obs = self.dm_5x5._subset_to_dataframe(["b", "d"], [])
        pdt.assert_frame_equal(obs, exp, check_dtype=False)
        obs = self.dm_5x5._subset_to_dataframe([], [])
        pdt.assert_frame_equal(obs, exp, check_dtype=False)

    def test_init_from_dm(self):
        ids = ["foo", "bar", "baz"]

        # PairwiseMatrix -> PairwiseMatrix
        exp = self.matobj(self.dm_3x3_data, ids)
        obs = self.matobj(self.dm_3x3, ids)
        self.assertEqual(obs, exp)
        # Test that copy of data is not made.
        self.assertTrue(obs.data is self.dm_3x3.data)
        obs.data[0, 1] = 424242
        self.assertTrue(np.array_equal(obs.data, self.dm_3x3.data))

        # DistanceMatrix -> PairwiseMatrix
        exp = self.matobj(self.dm_3x3_data, ids)
        obs = self.matobj(self.matobj(self.dm_3x3_data, ("a", "b", "c")), ids)
        self.assertEqual(obs, exp)

        # PairwiseMatrix -> DistanceMatrix
        with self.assertRaises(DistanceMatrixError):
            DistanceMatrix(self.dm_2x2_asym, ["foo", "bar"])

    def test_init_non_hollow_dm(self):
        data = [[1, 1], [1, 1]]
        obs = self.matobj(data, ["a", "b"])
        self.assertTrue(np.array_equal(obs.data, data))
        data_hollow = skbio.stats.distance._utils.is_hollow(obs.data)
        self.assertEqual(data_hollow, False)

    def test_init_no_ids(self):
        exp = self.matobj(self.dm_3x3_data, ("0", "1", "2"))
        obs = self.matobj(self.dm_3x3_data)
        self.assertEqual(obs, exp)
        self.assertEqual(obs["1", "2"], 12.0)

    def test_init_invalid_input(self):
        # Empty data.
        with self.assertRaises(PairwiseMatrixError):
            self.matobj([], [])

        # Another type of empty data.
        with self.assertRaises(PairwiseMatrixError):
            self.matobj(np.empty((0, 0)), [])

        # Invalid number of dimensions.
        with self.assertRaises(PairwiseMatrixError):
            self.matobj([1, 2, 3], ["a"])

        # Dimensions don't match.
        with self.assertRaises(PairwiseMatrixError):
            self.matobj([[1, 2, 3]], ["a"])

        data = [[0, 1], [1, 0]]

        # Duplicate IDs.
        with self.assertRaises(PairwiseMatrixError):
            self.matobj(data, ["a", "a"])

        # Number of IDs don't match dimensions.
        with self.assertRaises(PairwiseMatrixError):
            self.matobj(data, ["a", "b", "c"])
        with self.assertRaises(PairwiseMatrixError):
            self.matobj(data, [])

    def test_from_iterable_non_hollow_data(self):
        iterable = (x for x in range(4))

        exp = self.matobj([[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]])
        res = self.matobj.from_iterable(iterable, lambda a, b: 1)
        self.assertEqual(res, exp)

    def test_from_iterable_asymmetric_data(self):
        iterable = (x for x in range(4))

        exp = self.matobj(
            [[0, 1, 2, 3], [-1, 0, 1, 2], [-2, -1, 0, 1], [-3, -2, -1, 0]]
        )
        res = self.matobj.from_iterable(iterable, lambda a, b: b - a)
        self.assertEqual(res, exp)

    def test_from_iterable_no_key(self):
        iterable = (x for x in range(4))

        exp = self.matobj([[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]])
        res = self.matobj.from_iterable(iterable, lambda a, b: abs(b - a))
        self.assertEqual(res, exp)

    def test_from_iterable_with_key(self):
        iterable = (x for x in range(4))

        exp = self.matobj(
            [[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]],
            ["0", "1", "4", "9"],
        )
        res = self.matobj.from_iterable(
            iterable, lambda a, b: abs(b - a), key=lambda x: str(x**2)
        )
        self.assertEqual(res, exp)

    def test_from_iterable_empty(self):
        with self.assertRaises(PairwiseMatrixError):
            self.matobj.from_iterable([], lambda x: x)

    def test_from_iterable_single(self):
        exp = self.matobj([[100]])
        res = self.matobj.from_iterable(["boo"], lambda a, b: 100)
        self.assertEqual(res, exp)

    def test_from_iterable_with_keys(self):
        iterable = (x for x in range(4))

        exp = self.matobj(
            [[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]],
            ["0", "1", "4", "9"],
        )
        res = self.matobj.from_iterable(
            iterable, lambda a, b: abs(b - a), keys=iter(["0", "1", "4", "9"])
        )
        self.assertEqual(res, exp)

    def test_from_iterable_with_key_and_keys(self):
        iterable = (x for x in range(4))
        with self.assertRaises(ValueError):
            self.matobj.from_iterable(
                iterable, lambda a, b: abs(b - a), key=str, keys=["1", "2", "3", "4"]
            )

    def test_from_iterable_scipy_hamming_metric_with_metadata(self):
        # test for #1254
        seqs = [
            Sequence("ACGT"),
            Sequence("ACGA", metadata={"id": "seq1"}),
            Sequence("AAAA", metadata={"id": "seq2"}),
            Sequence("AAAA", positional_metadata={"qual": range(4)}),
        ]

        exp = self.matobj(
            [
                [0, 0.25, 0.75, 0.75],
                [0.25, 0.0, 0.5, 0.5],
                [0.75, 0.5, 0.0, 0.0],
                [0.75, 0.5, 0.0, 0.0],
            ],
            ["a", "b", "c", "d"],
        )

        dm = self.matobj.from_iterable(
            seqs, metric=scipy.spatial.distance.hamming, keys=["a", "b", "c", "d"]
        )

        self.assertEqual(dm, exp)

    def test_from_iterable_skbio_hamming_metric_with_metadata(self):
        # test for #1254
        seqs = [
            Sequence("ACGT"),
            Sequence("ACGA", metadata={"id": "seq1"}),
            Sequence("AAAA", metadata={"id": "seq2"}),
            Sequence("AAAA", positional_metadata={"qual": range(4)}),
        ]

        exp = self.matobj(
            [
                [0, 0.25, 0.75, 0.75],
                [0.25, 0.0, 0.5, 0.5],
                [0.75, 0.5, 0.0, 0.0],
                [0.75, 0.5, 0.0, 0.0],
            ],
            ["a", "b", "c", "d"],
        )

        dm = self.matobj.from_iterable(
            seqs, metric=skbio.sequence.distance.hamming, keys=["a", "b", "c", "d"]
        )

        self.assertEqual(dm, exp)

    def test_data(self):
        for dm, exp in zip(self.dms, self.dm_redundant_forms):
            obs = dm.data
            self.assertTrue(np.array_equal(obs, exp))

        with self.assertRaises(AttributeError):
            self.dm_3x3.data = "foo"

    def test_ids(self):
        obs = self.dm_3x3.ids
        self.assertEqual(obs, ("a", "b", "c"))

        # Test that we overwrite the existing IDs and that the ID index is
        # correctly rebuilt.
        new_ids = ["foo", "bar", "baz"]
        self.dm_3x3.ids = new_ids
        obs = self.dm_3x3.ids
        self.assertEqual(obs, tuple(new_ids))
        self.assertTrue(np.array_equal(self.dm_3x3["bar"], np.array([0.01, 0.0, 12.0])))
        with self.assertRaises(MissingIDError):
            self.dm_3x3["b"]

    def test_ids_invalid_input(self):
        with self.assertRaises(PairwiseMatrixError):
            self.dm_3x3.ids = ["foo", "bar"]
        # Make sure that we can still use the dissimilarity matrix after trying
        # to be evil.
        obs = self.dm_3x3.ids
        self.assertEqual(obs, ("a", "b", "c"))

    def test_ids_wrong_size(self):
        with self.assertRaises(PairwiseMatrixError):
            PairwiseMatrix([1, 2, 3], ["a"])

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
        self.assertEqual(self.dm_3x3.index("a"), 0)
        self.assertEqual(self.dm_3x3.index("b"), 1)
        self.assertEqual(self.dm_3x3.index("c"), 2)

        with self.assertRaises(MissingIDError):
            self.dm_3x3.index("d")

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

        new_ids = ["hello", "world"]
        copy.ids = new_ids
        self.assertNotEqual(copy.ids, self.dm_2x2.ids)

        copy = self.dm_2x2.copy()
        copy.data[0, 1] = 0.0001
        self.assertFalse(np.array_equal(copy.data, self.dm_2x2.data))

    def test_filter_no_filtering(self):
        # Don't actually filter anything -- ensure we get back a different
        # object.
        obs = self.dm_3x3.filter(["a", "b", "c"])
        self.assertEqual(obs, self.dm_3x3)
        self.assertFalse(obs is self.dm_3x3)

    def test_filter_reorder(self):
        # Don't filter anything, but reorder the distance matrix.
        order = ["c", "a", "b"]
        exp = self.matobj([[0, 4.2, 12], [4.2, 0, 0.01], [12, 0.01, 0]], order)
        obs = self.dm_3x3.filter(order)
        self.assertEqual(obs, exp)

    def test_filter_single_id(self):
        ids = ["b"]
        exp = self.matobj([[0]], ids)
        obs = self.dm_2x2_asym.filter(ids)
        self.assertEqual(obs, exp)

    def test_filter_asymmetric(self):
        # 2x2
        ids = ["b", "a"]
        exp = self.matobj([[0, -2], [1, 0]], ids)
        obs = self.dm_2x2_asym.filter(ids)
        self.assertEqual(obs, exp)

        # 3x3
        dm = self.matobj(
            [[0, 10, 53], [42, 0, 22.5], [53, 1, 0]], ("bro", "brah", "breh")
        )
        ids = ["breh", "brah"]
        exp = self.matobj([[0, 1], [22.5, 0]], ids)
        obs = dm.filter(ids)
        self.assertEqual(obs, exp)

    def test_filter_subset(self):
        ids = ("c", "a")
        exp = self.matobj([[0, 4.2], [4.2, 0]], ids)
        obs = self.dm_3x3.filter(ids)
        self.assertEqual(obs, exp)

        ids = ("b", "a")
        exp = self.matobj([[0, 0.01], [0.01, 0]], ids)
        obs = self.dm_3x3.filter(ids)
        self.assertEqual(obs, exp)

        # 4x4
        dm = self.matobj([[0, 1, 55, 7], [1, 0, 16, 1], [55, 16, 0, 23], [7, 1, 23, 0]])
        ids = np.asarray(["3", "0", "1"])
        exp = self.matobj([[0, 7, 1], [7, 0, 1], [1, 1, 0]], ids)
        obs = dm.filter(ids)
        self.assertEqual(obs, exp)

    def test_filter_duplicate_ids(self):
        with self.assertRaises(PairwiseMatrixError):
            self.dm_3x3.filter(["c", "a", "c"])

    def test_filter_missing_ids(self):
        with self.assertRaises(MissingIDError):
            self.dm_3x3.filter(["c", "bro"])

    def test_filter_missing_ids_strict_false(self):
        # no exception should be raised
        ids = ("c", "a")
        exp = self.matobj([[0, 4.2], [4.2, 0]], ids)
        obs = self.dm_3x3.filter(["c", "a", "not found"], strict=False)
        self.assertEqual(obs, exp)

    def test_filter_empty_ids(self):
        with self.assertRaises(PairwiseMatrixError):
            self.dm_3x3.filter([])

    def test_filter_sparse_same_ids(self):
        pm = self.matobj([1, 0, 0], ids=["a", "b", "c"], sparse=True)
        pmf = pm.filter(["a", "b", "c"])
        self.assertTrue(pm._flags["SPARSE"] and pmf._flags["SPARSE"])
        self.assertTrue(issparse(pm.data) and issparse(pmf.data))
        npt.assert_array_equal(pm.data.todense(), pmf.data.todense())


    @skipUnless(has_matplotlib, "Matplotlib not available.")
    def test_plot_default(self):
        fig = self.dm_1x1.plot()
        self.assertIsInstance(fig, mpl.figure.Figure)
        axes = fig.get_axes()
        self.assertEqual(len(axes), 2)
        ax = axes[0]
        self.assertEqual(ax.get_title(), "")
        xticks = []
        for tick in ax.get_xticklabels():
            xticks.append(tick.get_text())
        self.assertEqual(xticks, ["a"])
        yticks = []
        for tick in ax.get_yticklabels():
            yticks.append(tick.get_text())
        self.assertEqual(yticks, ["a"])

    @skipUnless(has_matplotlib, "Matplotlib not available.")
    def test_plot_no_default(self):
        ids = ["0", "one", "2", "three", "4.000"]
        data = (
            [0, 1, 2, 3, 4],
            [1, 0, 1, 2, 3],
            [2, 1, 0, 1, 2],
            [3, 2, 1, 0, 1],
            [4, 3, 2, 1, 0],
        )
        dm = self.matobj(data, ids)
        fig = dm.plot(cmap="Reds", title="Testplot")
        self.assertIsInstance(fig, mpl.figure.Figure)
        axes = fig.get_axes()
        self.assertEqual(len(axes), 2)
        ax = axes[0]
        self.assertEqual(ax.get_title(), "Testplot")
        xticks = []
        for tick in ax.get_xticklabels():
            xticks.append(tick.get_text())
        self.assertEqual(xticks, ["0", "one", "2", "three", "4.000"])
        yticks = []
        for tick in ax.get_yticklabels():
            yticks.append(tick.get_text())
        self.assertEqual(yticks, ["0", "one", "2", "three", "4.000"])

    def test_to_data_frame_1x1(self):
        df = self.dm_1x1.to_data_frame()
        exp = pd.DataFrame([[0.0]], index=["a"], columns=["a"])
        assert_data_frame_almost_equal(df, exp)

    def test_to_data_frame_3x3(self):
        df = self.dm_3x3.to_data_frame()
        exp = pd.DataFrame(
            [[0.0, 0.01, 4.2], [0.01, 0.0, 12.0], [4.2, 12.0, 0.0]],
            index=["a", "b", "c"],
            columns=["a", "b", "c"],
        )
        assert_data_frame_almost_equal(df, exp)

    def test_to_data_frame_default_ids(self):
        df = self.matobj(self.dm_2x2_data).to_data_frame()
        exp = pd.DataFrame(
            [[0.0, 0.123], [0.123, 0.0]], index=["0", "1"], columns=["0", "1"]
        )
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
        self.assertTrue(self.dm_3x3 != "foo")

        # Wrong shape.
        self.assertTrue(self.dm_3x3 != self.dm_1x1)

        # Wrong IDs.
        other = self.dm_3x3.copy()
        other.ids = ["foo", "bar", "baz"]
        self.assertTrue(self.dm_3x3 != other)

        # Wrong data.
        other = self.dm_3x3.copy()
        other.data[1, 0] = 42.42
        self.assertTrue(self.dm_3x3 != other)

        self.assertFalse(self.dm_2x2 != self.dm_2x2)

    def test_contains(self):
        self.assertTrue("a" in self.dm_3x3)
        self.assertTrue("b" in self.dm_3x3)
        self.assertTrue("c" in self.dm_3x3)
        self.assertFalse("d" in self.dm_3x3)

    def test_getslice(self):
        # Slice of first dimension only. Test that __getslice__ defers to
        # __getitem__.
        obs = self.dm_2x2[1:]
        self.assertTrue(np.array_equal(obs, np.array([[0.123, 0.0]])))
        self.assertEqual(type(obs), np.ndarray)

    def test_getitem_by_id(self):
        obs = self.dm_1x1["a"]
        self.assertTrue(np.array_equal(obs, np.array([0.0])))

        obs = self.dm_2x2_asym["b"]
        self.assertTrue(np.array_equal(obs, np.array([-2.0, 0.0])))

        obs = self.dm_3x3["c"]
        self.assertTrue(np.array_equal(obs, np.array([4.2, 12.0, 0.0])))

        with self.assertRaises(MissingIDError):
            self.dm_2x2["c"]

    def test_getitem_by_id_pair(self):
        # Same object.
        self.assertEqual(self.dm_1x1["a", "a"], 0.0)

        # Different objects (symmetric).
        self.assertEqual(self.dm_3x3["b", "c"], 12.0)
        self.assertEqual(self.dm_3x3["c", "b"], 12.0)

        # Different objects (asymmetric).
        self.assertEqual(self.dm_2x2_asym["a", "b"], 1.0)
        self.assertEqual(self.dm_2x2_asym["b", "a"], -2.0)

        with self.assertRaises(MissingIDError):
            self.dm_2x2["a", "c"]

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
        with self.assertRaises(PairwiseMatrixError):
            self.dm_3x3._validate_shape(np.array([[0, 42], [42, 0]]))

    def test_validate_invalid_shape(self):
        # first check it actually likes good matrices
        self.dm_3x3._validate_shape(np.array([[0.0, 42.0], [42.0, 0.0]]))
        # it checks just the shape, not the content
        self.dm_3x3._validate_shape(np.array([[1.0, 2.0], [3.0, 4.0]]))
        # empty array
        with self.assertRaises(PairwiseMatrixError):
            self.dm_3x3._validate_shape(np.array([]))
        # invalid shape
        with self.assertRaises(PairwiseMatrixError):
            self.dm_3x3._validate_shape(
                np.array([[0.0, 42.0], [42.0, 0.0], [22.0, 22.0]])
            )
        with self.assertRaises(PairwiseMatrixError):
            self.dm_3x3._validate_shape(
                np.array([[[0.0, 42.0], [42.0, 0.0]], [[0.0, 24.0], [24.0, 0.0]]])
            )

    def test_validate_invalid_shape_sparse(self):
        # first check it actually likes good matrices
        self.dm_3x3._validate_shape(csr_array(np.array([[0.0, 42.0], [42.0, 0.0]])))
        # it checks just the shape, not the content
        self.dm_3x3._validate_shape(csr_array(np.array([[1.0, 2.0], [3.0, 4.0]])))
        # empty array
        with self.assertRaises(PairwiseMatrixError):
            self.dm_3x3._validate_shape(csr_array(np.array([])))
        # invalid shape
        with self.assertRaises(PairwiseMatrixError):
            self.dm_3x3._validate_shape(
                csr_array(np.array([[0.0, 42.0], [42.0, 0.0], [22.0, 22.0]]))
            )
        with self.assertRaises(PairwiseMatrixError):
            self.dm_3x3._validate_shape(
                csr_array(np.array([[0, 42], [42, 0], [22, 22]]))
            )

    def test_validate_invalid_ids(self):
        # repeated ids
        with self.assertRaises(PairwiseMatrixError):
            self.dm_3x3._validate_ids(self.dm_3x3.data, ["a", "a"])
        # empty ids
        with self.assertRaises(PairwiseMatrixError):
            self.dm_3x3._validate_ids(self.dm_3x3.data, [])
        # invalid shape
        with self.assertRaises(PairwiseMatrixError):
            self.dm_3x3._validate_ids(self.dm_3x3.data, ["a", "b", "c", "d"])

    def test_init__data_with_sparse(self):
        mat = csr_array([[0, 0, 1], [0, 1, 0], [0, 0, 0]], dtype=np.float64)
        pm = self.matobj(mat, sparse=True)
        self.assertTrue(pm._flags["SPARSE"])
        self.assertTrue(issparse(pm.data))
        npt.assert_array_equal(mat.todense(), pm.data.todense())

    def test_init_data_condensed_to_sparse(self):
        mat = np.array([1, 2, 3, 4, 5, 6])
        pm = self.matobj(mat, sparse=True)
        self.assertTrue(pm._flags["SPARSE"])
        self.assertTrue(issparse(pm.data))
        npt.assert_array_equal(pm.data.todense(), np.array([[0.0, 1.0, 2.0, 3.0],
                                                            [1.0, 0.0, 4.0, 5.0],
                                                            [2.0, 4.0, 0.0, 6.0],
                                                            [3.0, 5.0, 6.0, 0.0]]))

    def test_init_data_np_to_sparse(self):
        mat = np.array([[0.0, 1.0, 0.0],
                        [0.5, 2.0, 0.0],
                        [0.0, 0.0, 0.0]])
        pm = self.matobj(mat, sparse=True)
        self.assertTrue(pm._flags["SPARSE"])
        self.assertTrue(issparse(pm.data))
        npt.assert_array_equal(pm.data.todense(), mat)

    def test_init_data_sparse_to_dense(self):
        mat = csr_array([[0, 3, 1],
                         [2, 5, 0],
                         [0, 0, 0]], dtype=float)
        pm = self.matobj(mat, sparse=False)
        self.assertTrue(not pm._flags["SPARSE"])
        self.assertTrue(not issparse(pm.data))
        npt.assert_array_equal(mat.todense(), pm.data)

class SymmetricMatrixTestBase(PairwiseMatrixTestData):
    @classmethod
    def get_matrix_class(cls):
        return None

    def setUp(self):
        super(SymmetricMatrixTestBase, self).setUp()
        self.matobj = self.get_matrix_class()
        self.dm_1x1 = self.matobj(self.dm_1x1_data, ["a"])
        self.dm_2x2 = self.matobj(self.dm_2x2_data, ["a", "b"])
        self.dm_3x3 = self.matobj(self.dm_3x3_data, ["a", "b", "c"])
        self.dm_1x1_cond = self.matobj(self.dm_1x1_data, ["a"], condensed=True)
        self.dm_2x2_cond = self.matobj(self.dm_2x2_data, ["a", "b"], condensed=True)
        self.dm_3x3_cond = self.matobj(
            self.dm_3x3_data, ["a", "b", "c"], condensed=True
        )
        self.dm_5x5 = self.matobj([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], list("abcde"))
        self.dm_5x5_cond = self.matobj(
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], list("abcde"), condensed=True
        )
        self.dm_5x5_cond_diag = self.matobj(
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            list("abcde"),
            condensed=True,
            diagonal=[90, 80, 70, 60, 50],
        )

        self.dms = [
            self.dm_1x1,
            self.dm_2x2,
            self.dm_3x3,
            self.dm_1x1_cond,
            self.dm_2x2_cond,
            self.dm_3x3_cond,
        ]
        self.dm_condensed_forms = [
            np.array([]),
            np.array([0.123]),
            np.array([0.01, 4.2, 12.0]),
        ] * 2

    def test_get_element_from_condensed(self):
        # check that it gets diagonal values correct
        exp = 70
        obs = self.dm_5x5_cond_diag["c", "c"]
        self.assertEqual(obs, exp)

        # check non diagonal values
        exp = 1
        obs = self.dm_5x5_cond_diag["a", "b"]

    def test_get_row_from_condensed(self):
        exp = np.array([90, 1, 2, 3, 4])
        obs = self.dm_5x5_cond_diag["a"]
        npt.assert_equal(obs, exp)

        exp = np.array([1, 80, 5, 6, 7])
        obs = self.dm_5x5_cond_diag["b"]
        npt.assert_equal(obs, exp)

        exp = np.array([4, 7, 9, 10, 50])
        obs = self.dm_5x5_cond_diag["e"]
        npt.assert_equal(obs, exp)

    def test_matrix_from_matrix(self):
        # symmetric from symmetric
        sm = SymmetricMatrix(self.dm_5x5)
        npt.assert_equal(sm.data, self.dm_5x5.data)
        self.assertEqual(sm.ids, self.dm_5x5.ids)
        npt.assert_equal(sm.diagonal, self.dm_5x5.diagonal)

        # symmetric from pairwise
        pm = PairwiseMatrix(self.dm_3x3)
        sm = self.matobj(pm)
        npt.assert_equal(sm.data, self.dm_3x3.data)
        self.assertEqual(sm.ids, self.dm_3x3.ids)
        npt.assert_equal(sm.diagonal, self.dm_3x3.diagonal)

        # symmetric from distance
        dm = DistanceMatrix(self.dm_3x3)
        sm = self.matobj(dm)
        npt.assert_equal(sm.data, self.dm_3x3.data)
        self.assertEqual(sm.ids, self.dm_3x3.ids)
        npt.assert_equal(sm.diagonal, self.dm_3x3.diagonal)

    def test_construction_sparse_and_condensed_raises_error(self):
        """Test that sparse=True and condensed=True raises error."""
        with self.assertRaisesRegex(
            DistanceMatrixError, "Cannot create sparse matrix in condensed form"
        ):
            DistanceMatrix(
                self.dm_3x3_data, sparse=True, condensed=True
            )

    def test_init_data_sparse_to_sparse(self):
        mat = csr_array([[0, 1], [1, 0]], dtype=float)
        sm = self.matobj(mat, sparse=False)
        self.assertTrue(not sm._flags["SPARSE"])
        self.assertTrue(not issparse(sm.data))
        npt.assert_array_equal(mat.todense(), sm.data)

    def test_init_data_condensed_to_sparse(self):
        mat = np.array([1, 2, 3], dtype=float)
        sm = self.matobj(mat, sparse=True)
        self.assertTrue(sm._flags["SPARSE"])
        self.assertTrue(issparse(sm.data))
        npt.assert_array_equal(sm.data.todense(), np.array([[0.0, 1.0, 2.0],
                                                            [1.0, 0.0, 3.0],
                                                            [2.0, 3.0, 0.0]]))

    def test_matrix_from_matrix_condensed(self):
        # symmetric from symmetric
        sm = SymmetricMatrix(self.dm_5x5_cond, condensed=True)
        npt.assert_equal(sm.data, self.dm_5x5_cond.data)
        self.assertEqual(sm.ids, self.dm_5x5_cond.ids)
        npt.assert_equal(sm.diagonal, self.dm_5x5_cond.diagonal)

        # symmetric from symmetric, change ids
        sm = SymmetricMatrix(
            self.dm_5x5_cond, ids=['aa', 'bb', 'cc', 'dd', 'ee'], condensed=True
        )
        npt.assert_equal(sm.data, self.dm_5x5_cond.data)
        self.assertEqual(sm.ids, ('aa', 'bb', 'cc', 'dd', 'ee'))
        npt.assert_equal(sm.diagonal, self.dm_5x5_cond.diagonal)

        # symmetric from pairwise
        pm = PairwiseMatrix(self.dm_3x3_cond)
        sm = self.matobj(pm, condensed=True)
        npt.assert_equal(sm.data, self.dm_3x3_cond.data)
        self.assertEqual(sm.ids, self.dm_3x3_cond.ids)
        npt.assert_equal(sm.diagonal, self.dm_3x3_cond.diagonal)

        # symmetric from pairwise, change ids
        pm = PairwiseMatrix(self.dm_3x3_cond)
        sm = self.matobj(pm, ids=['aa', 'bb', 'cc'], condensed=True)
        npt.assert_equal(sm.data, self.dm_3x3_cond.data)
        self.assertEqual(sm.ids, ('aa', 'bb', 'cc'))
        npt.assert_equal(sm.diagonal, self.dm_3x3_cond.diagonal)

        # symmetric from distance
        dm = DistanceMatrix(self.dm_3x3_cond)
        sm = self.matobj(dm, condensed=True)
        npt.assert_equal(sm.data, self.dm_3x3_cond.data)
        self.assertEqual(sm.ids, self.dm_3x3_cond.ids)
        npt.assert_equal(sm.diagonal, self.dm_3x3_cond.diagonal)

        # symmetric from distance, change ids
        dm = DistanceMatrix(self.dm_3x3_cond)
        sm = self.matobj(dm, ids=['aa', 'bb', 'cc'], condensed=True)
        npt.assert_equal(sm.data, self.dm_3x3_cond.data)
        self.assertEqual(sm.ids, ('aa', 'bb', 'cc'))
        npt.assert_equal(sm.diagonal, self.dm_3x3_cond.diagonal)

        # symmetric from distance, change diagonal, only possible when the
        # original matrix is condensed, because we can't provide diagonal with a 2D
        # matrix as input
        dm = DistanceMatrix(self.dm_3x3_cond, condensed=True)
        sm = self.matobj(dm, dm.ids, diagonal=[11, 22, 33], condensed=True)
        npt.assert_equal(sm.data, self.dm_3x3_cond.data)
        self.assertEqual(sm.ids, self.dm_3x3_cond.ids)
        npt.assert_equal(sm.diagonal, np.array([11, 22, 33]))

        #
        dm = SymmetricMatrix(self.dm_3x3_cond, condensed=True)
        sm = self.matobj(dm, dm.ids, diagonal=[11, 22, 33], condensed=True)
        npt.assert_equal(sm.data, self.dm_3x3_cond.data)
        self.assertEqual(sm.ids, self.dm_3x3_cond.ids)
        npt.assert_equal(sm.diagonal, np.array([11, 22, 33]))

        # this should raise error because can't pass diagonal with 2D input
        with self.assertRaises(SymmetricMatrixError) as e:
            pm = PairwiseMatrix(self.dm_3x3_cond)
            sm = self.matobj(pm, pm.ids, diagonal=[11, 22, 33], condensed=True)
        self.assertEqual(
            str(e.exception),
            "Cannot provide diagonal when data matrix is 2D. Information "
            "contained along diagonal is ambiguous.",
        )

    def test_as_condensed(self):
        sm = self.dm_3x3.as_condensed()
        # data should be the same, but copied
        self.assertEqual(sm, self.dm_3x3_cond)
        self.assertFalse(sm.data is self.dm_3x3.data)
        # ids aren't copied because they are immutable
        self.assertTrue(sm.ids is self.dm_3x3.ids)
        self.assertFalse(sm.diagonal is self.dm_3x3.diagonal)

    def test_as_redundant(self):
        sm_ = self.dm_3x3_cond.as_redundant()
        self.assertEqual(sm_, self.dm_3x3)
        self.assertFalse(sm_.data is self.dm_3x3_cond.data)
        # ids aren't copied because they are immutable
        self.assertTrue(sm_.ids is self.dm_3x3_cond.ids)
        self.assertFalse(sm_.diagonal is self.dm_3x3_cond.diagonal)

    def test_validate_ids_1d(self):
        with self.assertRaises(PairwiseMatrixError) as e:
            SymmetricMatrix([1, 2, 3], ["a"], condensed=True)
        self.assertEqual(
            str(e.exception),
            "The number of IDs (1) must match the number of rows/columns in "
            "the data (3).",
        )

    def test_init_diagonal_1d_none(self):
        # diagonal=None and 1d input should return None, because the
        # underlying data structure is redundant
        sm = SymmetricMatrix([1, 2, 3])
        obs = sm.diagonal
        exp = None
        self.assertEqual(obs, exp)

    def test_init_diagonal_2d_zero_trace(self):
        # diagonal=None and 2d input should return None because it is redundant
        sm = SymmetricMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
        obs = sm.diagonal
        exp = None
        self.assertEqual(obs, exp)

    def test_init_diagonal_2d_nonzero_trace(self):
        # diagonal=None and 2d input should return None if it contains
        # non-zeros, because the underlying data structure is redundant
        sm = SymmetricMatrix([[1, 2, 3], [2, 3, 4], [3, 4, 5]])
        obs = sm.diagonal
        exp = None
        self.assertEqual(obs, exp)

    def test_init_diagonal_condensed(self):
        sm = SymmetricMatrix([[1, 2, 3], [2, 3, 4], [3, 4, 5]], condensed=True)
        exp = np.array([1, 3, 5])
        obs = sm.diagonal
        npt.assert_equal(obs, exp)

    def test_init_diagonal_given_scalar(self):
        sm = SymmetricMatrix([1, 2, 3], diagonal=1.5)
        obs = sm.diagonal
        exp = 1.5
        self.assertEqual(obs, exp)

    def test_init_diagonal_given_arraylike(self):
        sm = SymmetricMatrix([1, 2, 3], diagonal=[2, 4, 5])
        obs = sm.diagonal
        exp = np.array([2, 4, 5])
        npt.assert_equal(obs, exp)

    def test_validate_data(self):
        with self.assertRaises(SymmetricMatrixError) as e:
            SymmetricMatrix([[1, 2], [4, 1]])
        self.assertEqual(
            str(e.exception), "Data must be symmetric and cannot contain NaNs."
        )

    def test_validate_diagonal_multidim(self):
        with self.assertRaises(SymmetricMatrixError) as e:
            SymmetricMatrix([1, 2, 3], diagonal=[[1, 2]])
        self.assertEqual(
            str(e.exception),
            "Diagonal must be 1 dimensional if it is an array. Found 2 dimensions.",
        )

    # def test_validate_diagonal_mismatch_2D(self):
    #     with self.assertRaises(SymmetricMatrixError) as e:
    #         SymmetricMatrix([[0, 1], [1, 0]], diagonal=[1, 2, 3])
    #     self.assertEqual(
    #         str(e.exception),
    #         "Length of diagonal (3) does not match the shape of the matrix (2, 2).",
    #     )

    def test_validate_diagonal_mismatch_1D(self):
        with self.assertRaises(SymmetricMatrixError) as e:
            SymmetricMatrix([1], diagonal=[1, 2, 3])
        self.assertEqual(
            str(e.exception),
            "Length of diagonal (3) does not match the shape of the matrix (2, 2).",
        )

    def test_validate_diagonal_providing_redundant(self):
        with self.assertRaises(SymmetricMatrixError) as e:
            sm = SymmetricMatrix([[0, 1], [1, 0]], diagonal=[2, 3])
        self.assertEqual(
            str(e.exception),
            "Cannot provide diagonal when data matrix is 2D. Information "
            "contained along diagonal is ambiguous.",
        )

    def test_validate_shape_invalid_1D_size(self):
        with self.assertRaises(SymmetricMatrixError) as e:
            SymmetricMatrix([1, 2, 3, 4])
        self.assertEqual(
            str(e.exception),
            "Incompatible vector size. It must be a binomial coefficient n "
            "choose 2 for some integer n >= 2.",
        )

    def test_validate_shape_3D(self):
        with self.assertRaises(SymmetricMatrixError) as e:
            SymmetricMatrix([[[1]]])
        self.assertEqual(
            str(e.exception),
            "Data must be have either 1 or 2 dimensions. Found 3 dimensions.",
        )

    def test_condensed_form_condensed(self):
        obs = SymmetricMatrix([1, 2, 3], condensed=True).condensed_form()
        exp = np.array([1, 2, 3])
        npt.assert_equal(obs, exp)

    def test_condensed_form_redundant(self):
        obs = SymmetricMatrix([1, 2, 3], condensed=False).condensed_form()
        exp = np.array([1, 2, 3])
        npt.assert_equal(obs, exp)

    def test_subset_to_dataframe(self):
        exp = pd.DataFrame(
            [
                ["b", "a", 1.0],
                ["b", "d", 6.0],
                ["b", "e", 7.0],
                ["d", "a", 3.0],
                ["d", "d", 0.0],
                ["d", "e", 10.0],
            ],
            columns=["i", "j", "value"],
        )

        obs = self.dm_5x5._subset_to_dataframe(["b", "d"], ["a", "d", "e"])
        pdt.assert_frame_equal(obs, exp)

        # and the empty edge cases
        exp = pd.DataFrame(
            [], columns=["i", "j", "value"], index=pd.RangeIndex(start=0, stop=0)
        )

        obs = self.dm_5x5._subset_to_dataframe([], ["a", "d", "e"])
        pdt.assert_frame_equal(obs, exp, check_dtype=False)
        obs = self.dm_5x5._subset_to_dataframe(["b", "d"], [])
        pdt.assert_frame_equal(obs, exp, check_dtype=False)
        obs = self.dm_5x5._subset_to_dataframe([], [])
        pdt.assert_frame_equal(obs, exp, check_dtype=False)

    def test_subset_to_dataframe_condensed(self):
        exp = pd.DataFrame(
            [
                ["b", "a", 1.0],
                ["b", "d", 6.0],
                ["b", "e", 7.0],
                ["d", "a", 3.0],
                ["d", "d", 0.0],
                ["d", "e", 10.0],
            ],
            columns=["i", "j", "value"],
        )

        obs = self.dm_5x5_cond._subset_to_dataframe(["b", "d"], ["a", "d", "e"])
        pdt.assert_frame_equal(obs, exp)

        # and the empty edge cases
        exp = pd.DataFrame(
            [], columns=["i", "j", "value"], index=pd.RangeIndex(start=0, stop=0)
        )

        obs = self.dm_5x5_cond._subset_to_dataframe([], ["a", "d", "e"])
        pdt.assert_frame_equal(obs, exp, check_dtype=False)
        obs = self.dm_5x5_cond._subset_to_dataframe(["b", "d"], [])
        pdt.assert_frame_equal(obs, exp, check_dtype=False)
        obs = self.dm_5x5_cond._subset_to_dataframe([], [])
        pdt.assert_frame_equal(obs, exp, check_dtype=False)

    def test_subset_to_dataframe_condensed_diagonal(self):
        exp = pd.DataFrame(
            [
                ["b", "a", 1.0],
                ["b", "d", 6.0],
                ["b", "e", 7.0],
                ["d", "a", 3.0],
                ["d", "d", 60.0],
                ["d", "e", 10.0],
            ],
            columns=["i", "j", "value"],
        )

        obs = self.dm_5x5_cond_diag._subset_to_dataframe(["b", "d"], ["a", "d", "e"])

        pdt.assert_frame_equal(obs, exp)

        # and the empty edge cases
        exp = pd.DataFrame(
            [], columns=["i", "j", "value"], index=pd.RangeIndex(start=0, stop=0)
        )

        obs = self.dm_5x5_cond_diag._subset_to_dataframe([], ["a", "d", "e"])
        pdt.assert_frame_equal(obs, exp, check_dtype=False)
        obs = self.dm_5x5_cond_diag._subset_to_dataframe(["b", "d"], [])
        pdt.assert_frame_equal(obs, exp, check_dtype=False)
        obs = self.dm_5x5_cond_diag._subset_to_dataframe([], [])
        pdt.assert_frame_equal(obs, exp, check_dtype=False)

    # test filtering on condensed forms
    def test_filter_no_filtering(self):
        # Don't actually filter anything -- ensure we get back a different
        # object.
        obs = self.dm_3x3_cond.filter(["a", "b", "c"])
        self.assertEqual(obs, self.dm_3x3_cond)
        self.assertFalse(obs is self.dm_3x3_cond)

    def test_filter_reorder(self):
        # Don't filter anything, but reorder the distance matrix.
        order = ["c", "a", "b"]
        exp = self.matobj(
            [[0, 4.2, 12], [4.2, 0, 0.01], [12, 0.01, 0]], order, condensed=True
        )
        obs = self.dm_3x3_cond.filter(order)
        self.assertEqual(obs, exp)

    def test_filter_single_id(self):
        ids = ["b"]
        exp = self.matobj([[0]], ids, condensed=True)
        obs = self.dm_2x2_cond.filter(ids)
        self.assertEqual(obs, exp)

    def test_filter_subset(self):
        ids = ("c", "a")
        exp = self.matobj([[0, 4.2], [4.2, 0]], ids, condensed=True)
        obs = self.dm_3x3_cond.filter(ids)
        self.assertEqual(obs, exp)

        ids = ("b", "a")
        exp = self.matobj([[0, 0.01], [0.01, 0]], ids, condensed=True)
        obs = self.dm_3x3_cond.filter(ids)
        self.assertEqual(obs, exp)

        # 4x4
        dm = self.matobj(
            [[0, 1, 55, 7], [1, 0, 16, 1], [55, 16, 0, 23], [7, 1, 23, 0]],
            condensed=True,
        )
        ids = np.asarray(["3", "0", "1"])
        exp = self.matobj([[0, 7, 1], [7, 0, 1], [1, 1, 0]], ids, condensed=True)
        obs = dm.filter(ids)
        self.assertEqual(obs, exp)

    def test_filter_duplicate_ids(self):
        with self.assertRaises(PairwiseMatrixError) as e:
            self.dm_3x3_cond.filter(["c", "a", "c"])
        self.assertEqual(
            str(e.exception),
            "IDs must be unique. Found the following duplicate IDs: 'c'",
        )

    def test_filter_missing_ids(self):
        with self.assertRaises(MissingIDError) as e:
            self.dm_3x3_cond.filter(["c", "bro"])
        self.assertEqual(str(e.exception), "The ID 'bro' is not in the matrix.")

    def test_filter_missing_ids_strict_false(self):
        # no exception should be raised
        ids = ("c", "a")
        exp = self.matobj([[0, 4.2], [4.2, 0]], ids, condensed=True)
        obs = self.dm_3x3_cond.filter(["c", "a", "not found"], strict=False)
        self.assertEqual(obs, exp)

    def test_filter_empty_ids(self):
        with self.assertRaises(PairwiseMatrixError) as e:
            self.dm_3x3_cond.filter([])
        self.assertEqual(str(e.exception), "IDs must be at least 1 in size.")

    def test_getslice(self):
        # Slice of first dimension only. Test that __getslice__ defers to
        # __getitem__.
        obs = self.dm_2x2_cond[1:]
        self.assertTrue(np.array_equal(obs, np.array([[0.123, 0.0]])))
        self.assertEqual(type(obs), np.ndarray)

    def test_getitem_by_id(self):
        obs = self.dm_1x1_cond["a"]
        self.assertTrue(np.array_equal(obs, np.array([0.0])))

        obs = self.dm_2x2_cond["b"]
        self.assertTrue(np.array_equal(obs, np.array([0.123, 0.0])))

        obs = self.dm_3x3_cond["c"]
        self.assertTrue(np.array_equal(obs, np.array([4.2, 12.0, 0.0])))

        with self.assertRaises(MissingIDError) as e:
            self.dm_2x2_cond["c"]
        self.assertEqual(str(e.exception), "The ID 'c' is not in the matrix.")

    def test_getitem_by_id_with_scalar_diag(self):
        # get row
        sm = self.matobj([1, 2, 3], ['a', 'b', 'c'], diagonal=8, condensed=True)
        obs = sm["a"]
        exp = np.array([8, 1, 2])
        npt.assert_equal(obs, exp)

        # get single value
        obs = sm["a", "a"]
        exp = 8.0
        npt.assert_equal(obs, exp)

    def test_getitem_by_id_pair(self):
        # Same object.
        self.assertEqual(self.dm_1x1_cond["a", "a"], 0.0)

        # Different objects (symmetric).
        self.assertEqual(self.dm_3x3_cond["b", "c"], 12.0)
        self.assertEqual(self.dm_3x3_cond["c", "b"], 12.0)

        # Different objects.
        self.assertEqual(self.dm_2x2_cond["a", "b"], 0.123)
        self.assertEqual(self.dm_2x2_cond["b", "a"], 0.123)

        with self.assertRaises(MissingIDError) as e:
            self.dm_2x2_cond["a", "c"]
        self.assertEqual(str(e.exception), "The ID 'c' is not in the matrix.")

    def test_getitem_ndarray_indexing(self):
        # Single element access.
        obs = self.dm_3x3_cond[0, 1]
        self.assertEqual(obs, 0.01)

        # Single element access (via two __getitem__ calls).
        obs = self.dm_3x3_cond[0][1]
        self.assertEqual(obs, 0.01)

        # Row access.
        obs = self.dm_3x3_cond[1]
        self.assertTrue(np.array_equal(obs, np.array([0.01, 0.0, 12.0])))
        self.assertEqual(type(obs), np.ndarray)

        # Grab all data.
        obs = self.dm_3x3_cond[:, :]
        self.assertTrue(np.array_equal(obs, self.dm_3x3.data))
        self.assertEqual(type(obs), np.ndarray)

        with self.assertRaises(IndexError) as e:
            self.dm_3x3_cond[:, 3]
        self.assertEqual(
            str(e.exception), "index 3 is out of bounds for axis 1 with size 3"
        )


class DistanceMatrixTestBase(PairwiseMatrixTestData):
    @classmethod
    def get_matrix_class(cls):
        return None

    def setUp(self):
        super(DistanceMatrixTestBase, self).setUp()
        self.matobj = self.get_matrix_class()
        self.dm_1x1 = self.matobj(self.dm_1x1_data, ["a"])
        self.dm_2x2 = self.matobj(self.dm_2x2_data, ["a", "b"])
        self.dm_3x3 = self.matobj(self.dm_3x3_data, ["a", "b", "c"])
        self.dm_1x1_cond = self.matobj(self.dm_1x1_data, ["a"], condensed=True)
        self.dm_2x2_cond = self.matobj(self.dm_2x2_data, ["a", "b"], condensed=True)
        self.dm_3x3_cond = self.matobj(
            self.dm_3x3_data, ["a", "b", "c"], condensed=True
        )

        self.dms = [
            self.dm_1x1,
            self.dm_2x2,
            self.dm_3x3,
            self.dm_1x1_cond,
            self.dm_2x2_cond,
            self.dm_3x3_cond,
        ]
        self.dm_condensed_forms = [
            np.array([]),
            np.array([0.123]),
            np.array([0.01, 4.2, 12.0]),
        ] * 2

    def test_matrix_from_matrix(self):
        # distance from distance
        sm = DistanceMatrix(self.dm_3x3)
        npt.assert_equal(sm.data, self.dm_3x3.data)
        self.assertEqual(sm.ids, self.dm_3x3.ids)
        npt.assert_equal(sm.diagonal, self.dm_3x3.diagonal)

        # distance from pairwise
        pm = PairwiseMatrix(self.dm_3x3)
        sm = self.matobj(pm)
        npt.assert_equal(sm.data, self.dm_3x3.data)
        self.assertEqual(sm.ids, self.dm_3x3.ids)
        npt.assert_equal(sm.diagonal, self.dm_3x3.diagonal)

        # distance from symmetric
        dm = SymmetricMatrix(self.dm_3x3)
        sm = self.matobj(dm)
        npt.assert_equal(sm.data, self.dm_3x3.data)
        self.assertEqual(sm.ids, self.dm_3x3.ids)
        npt.assert_equal(sm.diagonal, self.dm_3x3.diagonal)

    def test_matrix_from_matrix_condensed(self):
        # distance from distance
        sm = DistanceMatrix(self.dm_3x3_cond, condensed=True)
        npt.assert_equal(sm.data, self.dm_3x3_cond.data)
        self.assertEqual(sm.ids, self.dm_3x3_cond.ids)
        npt.assert_equal(sm.diagonal, self.dm_3x3_cond.diagonal)

        # distance from pairwise
        pm = PairwiseMatrix(self.dm_3x3_cond)
        sm = self.matobj(pm, condensed=True)
        npt.assert_equal(sm.data, self.dm_3x3_cond.data)
        self.assertEqual(sm.ids, self.dm_3x3_cond.ids)
        npt.assert_equal(sm.diagonal, self.dm_3x3_cond.diagonal)

        # distance from symmetric
        dm = SymmetricMatrix(self.dm_3x3_cond)
        sm = self.matobj(dm, condensed=True)
        npt.assert_equal(sm.data, self.dm_3x3_cond.data)
        self.assertEqual(sm.ids, self.dm_3x3_cond.ids)
        npt.assert_equal(sm.diagonal, self.dm_3x3_cond.diagonal)

    def test_init_from_condensed_form(self):
        data = [1, 2, 3]
        exp = self.matobj([[0, 1, 2], [1, 0, 3], [2, 3, 0]], ["0", "1", "2"])
        res = self.matobj(data)
        self.assertEqual(exp, res)

    def test_init_from_condensed_to_condensed(self):
        data = [1, 2, 3]
        exp = self.matobj(
            [[0, 1, 2], [1, 0, 3], [2, 3, 0]], ["0", "1", "2"], condensed=True
        )
        res = self.matobj(data, condensed=True)
        self.assertEqual(exp, res)

    def test_init_invalid_input(self):
        # Asymmetric.
        data = [[0.0, 2.0], [1.0, 0.0]]
        with self.assertRaises(DistanceMatrixError) as e:
            self.matobj(data, ["a", "b"])
        self.assertEqual(
            str(e.exception), "Data must be symmetric and cannot contain NaNs."
        )

        # Non-hollow
        data = [[1.0, 2.0], [2.0, 1.0]]
        with self.assertRaises(DistanceMatrixError) as e:
            self.matobj(data, ["a", "b"])
        self.assertEqual(
            str(e.exception),
            "Data must  be hollow (i.e., the diagonal can only contain zeros).",
        )

        # Ensure that the superclass validation is still being performed.
        with self.assertRaises(PairwiseMatrixError) as e:
            self.matobj([[1, 2, 3]], ["a"])
        self.assertEqual(
            str(e.exception),
            "Data must be square (i.e., have the same number of rows and columns).",
        )

    def test_init_invalid_input_condensed(self):
        # Asymmetric.
        data = [[0.0, 2.0], [1.0, 0.0]]
        with self.assertRaises(DistanceMatrixError) as e:
            self.matobj(data, ["a", "b"], condensed=True)
        self.assertEqual(
            str(e.exception), "Data must be symmetric and cannot contain NaNs."
        )

        # Non-hollow
        data = [[1.0, 2.0], [2.0, 1.0]]
        with self.assertRaises(DistanceMatrixError) as e:
            self.matobj(data, ["a", "b"], condensed=True)
        self.assertEqual(
            str(e.exception),
            "Data must  be hollow (i.e., the diagonal can only contain zeros).",
        )

        # Ensure that the superclass validation is still being performed.
        with self.assertRaises(PairwiseMatrixError) as e:
            self.matobj([[1, 2, 3]], ["a"], condensed=True)
        self.assertEqual(
            str(e.exception),
            "Data must be square (i.e., have the same number of rows and columns).",
        )

    def test_init_nans(self):
        with self.assertRaisesRegex(DistanceMatrixError, r"NaNs"):
            self.matobj([[0.0, np.nan], [np.nan, 0.0]], ["a", "b"])

    def test_init_nans_condensed(self):
        with self.assertRaisesRegex(DistanceMatrixError, r"NaNs"):
            self.matobj([[0.0, np.nan], [np.nan, 0.0]], ["a", "b"], condensed=True)

    def test_from_iterable_no_key(self):
        iterable = (x for x in range(4))

        exp = self.matobj([[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]])
        res = self.matobj.from_iterable(iterable, lambda a, b: abs(b - a))
        self.assertEqual(res, exp)

    def test_from_iterable_no_key_condensed(self):
        iterable = (x for x in range(4))

        exp = self.matobj(
            [[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]], condensed=True
        )
        res = self.matobj.from_iterable(
            iterable, lambda a, b: abs(b - a), condensed=True
        )
        self.assertEqual(res, exp)

    def test_from_iterable_validate_equal_valid_data(self):
        validate_true = self.matobj.from_iterable(
            (x for x in range(4)), lambda a, b: abs(b - a), validate=True
        )
        validate_false = self.matobj.from_iterable(
            (x for x in range(4)), lambda a, b: abs(b - a), validate=False
        )
        self.assertEqual(validate_true, validate_false)

    def test_from_iterable_validate_equal_valid_data_condensed(self):
        validate_true = self.matobj.from_iterable(
            (x for x in range(4)),
            lambda a, b: abs(b - a),
            validate=True,
            condensed=True,
        )
        validate_false = self.matobj.from_iterable(
            (x for x in range(4)),
            lambda a, b: abs(b - a),
            validate=False,
            condensed=True,
        )
        self.assertEqual(validate_true, validate_false)

    def test_from_iterable_validate_false(self):
        iterable = (x for x in range(4))

        exp = self.matobj([[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]])
        res = self.matobj.from_iterable(
            iterable, lambda a, b: abs(b - a), validate=False
        )
        self.assertEqual(res, exp)

    def test_from_iterable_validate_false_condensed(self):
        iterable = (x for x in range(4))

        exp = self.matobj(
            [[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]], condensed=True
        )
        res = self.matobj.from_iterable(
            iterable, lambda a, b: abs(b - a), validate=False, condensed=True
        )
        self.assertEqual(res, exp)

    def test_from_iterable_with_invalid_diagonal(self):
        iterable = (x for x in range(4))
        with self.assertRaises(DistanceMatrixError) as e:
            self.matobj.from_iterable(
                iterable, lambda a, b: a * b, diagonal=1, validate=False
            )
        self.assertEqual(
            str(e.exception),
            "Data must  be hollow (i.e., the diagonal can only contain zeros)."
        )

    def test_from_iterable_with_valid_diagonal(self):
        iterable = (x for x in range(4))
        exp = self.matobj(
            [[0, 1, 2, 3], [1, 0, 2, 3], [2, 2, 0, 3], [3, 3, 3, 0]]
        )
        res = self.matobj.from_iterable(
            iterable, lambda a, b: a, diagonal=0, validate=False
        )
        self.assertTrue((res.data == exp.data).all())

    def test_from_iterable_validate_non_hollow(self):
        iterable = (x for x in range(4))
        with self.assertRaises(DistanceMatrixError) as e:
            self.matobj.from_iterable(iterable, lambda a, b: 1)
        self.assertEqual(
            str(e.exception),
            "Data must  be hollow (i.e., the diagonal can only contain zeros).",
        )

    def test_from_iterable_validate_non_hollow_condensed(self):
        iterable = (x for x in range(4))
        with self.assertRaises(DistanceMatrixError) as e:
            self.matobj.from_iterable(iterable, lambda a, b: 1, condensed=True)
        self.assertEqual(
            str(e.exception),
            "Data must  be hollow (i.e., the diagonal can only contain zeros).",
        )

    def test_from_iterable_validate_false_non_symmetric(self):
        exp = self.matobj([[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]])
        res = self.matobj.from_iterable(
            (x for x in range(4)), lambda a, b: a - b, validate=False
        )
        self.assertEqual(res, exp)

    def test_from_iterable_validate_false_non_symmetric_condensed(self):
        exp = self.matobj(
            [[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]], condensed=True
        )
        res = self.matobj.from_iterable(
            (x for x in range(4)), lambda a, b: a - b, validate=False, condensed=True
        )
        self.assertEqual(res, exp)

    def test_from_iterable_validate_asym(self):
        iterable = (x for x in range(4))
        with self.assertRaises(DistanceMatrixError):
            self.matobj.from_iterable(iterable, lambda a, b: b - a)

    def test_from_iterable_validate_asym_condensed(self):
        iterable = (x for x in range(4))
        with self.assertRaises(DistanceMatrixError):
            self.matobj.from_iterable(iterable, lambda a, b: b - a, condensed=True)

    def test_from_iterable_with_key(self):
        iterable = (x for x in range(4))

        exp = self.matobj(
            [[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]],
            ["0", "1", "4", "9"],
        )
        res = self.matobj.from_iterable(
            iterable, lambda a, b: abs(b - a), key=lambda x: str(x**2)
        )
        self.assertEqual(res, exp)

    def test_from_iterable_with_key_condensed(self):
        iterable = (x for x in range(4))

        exp = self.matobj(
            [[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]],
            ["0", "1", "4", "9"],
            condensed=True,
        )
        res = self.matobj.from_iterable(
            iterable, lambda a, b: abs(b - a), key=lambda x: str(x**2), condensed=True
        )
        self.assertEqual(res, exp)

    def test_from_iterable_empty(self):
        with self.assertRaises(PairwiseMatrixError):
            self.matobj.from_iterable([], lambda x: x)

    def test_from_iterable_empty_condensed(self):
        with self.assertRaises(PairwiseMatrixError):
            self.matobj.from_iterable([], lambda x: x, condensed=True)

    def test_from_iterable_single(self):
        exp = self.matobj([[0]])
        res = self.matobj.from_iterable(["boo"], lambda a, b: 0)
        self.assertEqual(res, exp)

    def test_from_iterable_single_condensed(self):
        exp = self.matobj([[0]], condensed=True)
        res = self.matobj.from_iterable(["boo"], lambda a, b: 0, condensed=True)
        self.assertEqual(res, exp)

    def test_from_iterable_with_keys(self):
        iterable = (x for x in range(4))

        exp = self.matobj(
            [[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]],
            ["0", "1", "4", "9"],
        )
        res = self.matobj.from_iterable(
            iterable, lambda a, b: abs(b - a), keys=iter(["0", "1", "4", "9"])
        )
        self.assertEqual(res, exp)

    def test_from_iterable_with_keys_condensed(self):
        iterable = (x for x in range(4))

        exp = self.matobj(
            [[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]],
            ["0", "1", "4", "9"],
            condensed=True,
        )
        res = self.matobj.from_iterable(
            iterable,
            lambda a, b: abs(b - a),
            keys=iter(["0", "1", "4", "9"]),
            condensed=True,
        )
        self.assertEqual(res, exp)

    def test_from_iterable_with_key_and_keys(self):
        iterable = (x for x in range(4))
        with self.assertRaises(ValueError):
            self.matobj.from_iterable(
                iterable, lambda a, b: abs(b - a), key=str, keys=["1", "2", "3", "4"]
            )

        iterable = (x for x in range(4))
        with self.assertRaises(ValueError):
            self.matobj.from_iterable(
                iterable,
                lambda a, b: abs(b - a),
                key=str,
                keys=["1", "2", "3", "4"],
                condensed=True,
            )

    def test_from_iterable_scipy_hamming_metric_with_metadata(self):
        # test for #1254
        seqs = [
            Sequence("ACGT"),
            Sequence("ACGA", metadata={"id": "seq1"}),
            Sequence("AAAA", metadata={"id": "seq2"}),
            Sequence("AAAA", positional_metadata={"qual": range(4)}),
        ]

        exp = self.matobj(
            [
                [0, 0.25, 0.75, 0.75],
                [0.25, 0.0, 0.5, 0.5],
                [0.75, 0.5, 0.0, 0.0],
                [0.75, 0.5, 0.0, 0.0],
            ],
            ["a", "b", "c", "d"],
        )

        dm = self.matobj.from_iterable(
            seqs, metric=scipy.spatial.distance.hamming, keys=["a", "b", "c", "d"]
        )

        self.assertEqual(dm, exp)

    def test_from_iterable_scipy_hamming_metric_with_metadata_condensed(self):
        # test for #1254
        seqs = [
            Sequence("ACGT"),
            Sequence("ACGA", metadata={"id": "seq1"}),
            Sequence("AAAA", metadata={"id": "seq2"}),
            Sequence("AAAA", positional_metadata={"qual": range(4)}),
        ]

        exp = self.matobj(
            [
                [0, 0.25, 0.75, 0.75],
                [0.25, 0.0, 0.5, 0.5],
                [0.75, 0.5, 0.0, 0.0],
                [0.75, 0.5, 0.0, 0.0],
            ],
            ["a", "b", "c", "d"],
            condensed=True,
        )

        dm = self.matobj.from_iterable(
            seqs,
            metric=scipy.spatial.distance.hamming,
            keys=["a", "b", "c", "d"],
            condensed=True,
        )

        self.assertEqual(dm, exp)

    def test_from_iterable_skbio_hamming_metric_with_metadata(self):
        # test for #1254
        seqs = [
            Sequence("ACGT"),
            Sequence("ACGA", metadata={"id": "seq1"}),
            Sequence("AAAA", metadata={"id": "seq2"}),
            Sequence("AAAA", positional_metadata={"qual": range(4)}),
        ]

        exp = self.matobj(
            [
                [0, 0.25, 0.75, 0.75],
                [0.25, 0.0, 0.5, 0.5],
                [0.75, 0.5, 0.0, 0.0],
                [0.75, 0.5, 0.0, 0.0],
            ],
            ["a", "b", "c", "d"],
        )

        dm = self.matobj.from_iterable(
            seqs, metric=skbio.sequence.distance.hamming, keys=["a", "b", "c", "d"]
        )

        self.assertEqual(dm, exp)

    def test_from_iterable_skbio_hamming_metric_with_metadata_condensed(self):
        # test for #1254
        seqs = [
            Sequence("ACGT"),
            Sequence("ACGA", metadata={"id": "seq1"}),
            Sequence("AAAA", metadata={"id": "seq2"}),
            Sequence("AAAA", positional_metadata={"qual": range(4)}),
        ]

        exp = self.matobj(
            [
                [0, 0.25, 0.75, 0.75],
                [0.25, 0.0, 0.5, 0.5],
                [0.75, 0.5, 0.0, 0.0],
                [0.75, 0.5, 0.0, 0.0],
            ],
            ["a", "b", "c", "d"],
            condensed=True,
        )

        dm = self.matobj.from_iterable(
            seqs,
            metric=skbio.sequence.distance.hamming,
            keys=["a", "b", "c", "d"],
            condensed=True,
        )

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

        obs = self.dm_3x3.permute(condensed=True, seed=3)
        npt.assert_equal(obs, np.array([12.0, 4.2, 0.01]))

        obs = self.dm_3x3.permute(condensed=True, seed=2)
        npt.assert_equal(obs, np.array([4.2, 12.0, 0.01]))

        # Ensure dm hasn't changed after calling permute() on it a couple of
        # times.
        self.assertEqual(self.dm_3x3, dm_copy)

    def test_permute_condensed_on_condensed(self):
        # Can't really permute a 1x1 or 2x2...
        for _ in range(2):
            obs = self.dm_1x1_cond.permute(condensed=True)
            npt.assert_equal(obs, np.array([]))

        for _ in range(2):
            obs = self.dm_2x2_cond.permute(condensed=True)
            npt.assert_equal(obs, np.array([0.123]))

        dm_copy = self.dm_3x3_cond.copy()

        obs = self.dm_3x3_cond.permute(condensed=True, seed=3)
        npt.assert_equal(obs, np.array([12.0, 4.2, 0.01]))

        obs = self.dm_3x3_cond.permute(condensed=True, seed=2)
        npt.assert_equal(obs, np.array([4.2, 12.0, 0.01]))

        # Ensure dm hasn't changed after calling permute() on it a couple of
        # times.
        self.assertEqual(self.dm_3x3_cond, dm_copy)

    def test_permute_not_condensed(self):
        obs = self.dm_1x1.permute()
        self.assertEqual(obs, self.dm_1x1)
        self.assertFalse(obs is self.dm_1x1)

        obs = self.dm_2x2.permute()
        self.assertEqual(obs, self.dm_2x2)
        self.assertFalse(obs is self.dm_2x2)

        exp = self.matobj(
            [[0, 12, 4.2], [12, 0, 0.01], [4.2, 0.01, 0]], self.dm_3x3.ids
        )
        obs = self.dm_3x3.permute(seed=3)
        self.assertEqual(obs, exp)

        exp = self.matobj(
            [[0, 4.2, 12], [4.2, 0, 0.01], [12, 0.01, 0]], self.dm_3x3.ids
        )
        obs = self.dm_3x3.permute(seed=2)
        self.assertEqual(obs, exp)

    def test_permute_not_condensed_on_condensed(self):
        obs = self.dm_1x1_cond.permute()
        self.assertEqual(obs, self.dm_1x1_cond)
        self.assertFalse(obs is self.dm_1x1_cond)

        obs = self.dm_2x2_cond.permute()
        self.assertEqual(obs, self.dm_2x2_cond)
        self.assertFalse(obs is self.dm_2x2_cond)

        exp = self.matobj(
            [[0, 12, 4.2], [12, 0, 0.01], [4.2, 0.01, 0]],
            self.dm_3x3_cond.ids,
            condensed=True,
        )
        obs = self.dm_3x3_cond.permute(seed=3)
        self.assertEqual(obs, exp)

        exp = self.matobj(
            [[0, 4.2, 12], [4.2, 0, 0.01], [12, 0.01, 0]],
            self.dm_3x3_cond.ids,
            condensed=True,
        )
        obs = self.dm_3x3_cond.permute(seed=2)
        self.assertEqual(obs, exp)

    # TODO: need to test this with SymmetricMatrix vs DistanceMatrix for condensed form
    def test_eq(self):
        # Compare DistanceMatrix to PairwiseMatrix, where both have the
        # same data and IDs.
        eq_dm = PairwiseMatrix(self.dm_3x3_data, ["a", "b", "c"])
        self.assertTrue(self.dm_3x3 == eq_dm)
        self.assertTrue(eq_dm == self.dm_3x3)

    def test_to_series_1x1(self):
        series = self.dm_1x1.to_series()

        exp = pd.Series([], index=[], dtype="float64")
        assert_series_almost_equal(series, exp)

    def test_to_series_1x1_condensed(self):
        series = self.dm_1x1_cond.to_series()

        exp = pd.Series([], index=[], dtype="float64")
        assert_series_almost_equal(series, exp)

    def test_to_series_2x2(self):
        series = self.dm_2x2.to_series()

        exp = pd.Series([0.123], index=pd.Index([("a", "b")]))
        assert_series_almost_equal(series, exp)

    def test_to_series_2x2_condensed(self):
        series = self.dm_2x2_cond.to_series()

        exp = pd.Series([0.123], index=pd.Index([("a", "b")]))
        assert_series_almost_equal(series, exp)

    def test_to_series_4x4(self):
        dm = self.matobj(
            [
                [0.0, 0.2, 0.3, 0.4],
                [0.2, 0.0, 0.5, 0.6],
                [0.3, 0.5, 0.0, 0.7],
                [0.4, 0.6, 0.7, 0.0],
            ],
            ["a", "b", "c", "d"],
        )

        series = dm.to_series()

        exp = pd.Series(
            [0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
            index=pd.Index(
                [("a", "b"), ("a", "c"), ("a", "d"), ("b", "c"), ("b", "d"), ("c", "d")]
            ),
        )
        assert_series_almost_equal(series, exp)

    def test_to_series_4x4_condensed(self):
        dm = self.matobj(
            [
                [0.0, 0.2, 0.3, 0.4],
                [0.2, 0.0, 0.5, 0.6],
                [0.3, 0.5, 0.0, 0.7],
                [0.4, 0.6, 0.7, 0.0],
            ],
            ["a", "b", "c", "d"],
            condensed=True,
        )

        series = dm.to_series()

        exp = pd.Series(
            [0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
            index=pd.Index(
                [("a", "b"), ("a", "c"), ("a", "d"), ("b", "c"), ("b", "d"), ("c", "d")]
            ),
        )
        assert_series_almost_equal(series, exp)

    def test_to_series_default_ids(self):
        series = self.matobj(self.dm_2x2_data).to_series()

        exp = pd.Series([0.123], index=pd.Index([("0", "1")]))
        assert_series_almost_equal(series, exp)

    def test_to_series_default_ids_condensed(self):
        series = self.matobj(self.dm_2x2_data, condensed=True).to_series()

        exp = pd.Series([0.123], index=pd.Index([("0", "1")]))
        assert_series_almost_equal(series, exp)

    def test_validate_asym_shape(self):
        # first check it actually likes good matrices
        data_good = np.array([[0.0, 42.0], [42.0, 0.0]])
        data_sym, data_hollow = is_symmetric_and_hollow(data_good)
        self.assertEqual(data_sym, True)
        del data_sym
        self.assertEqual(data_hollow, True)
        del data_hollow
        data_sym = skbio.stats.distance._utils.is_symmetric(data_good)
        self.assertEqual(data_sym, True)
        del data_sym
        data_hollow = skbio.stats.distance._utils.is_hollow(data_good)
        self.assertEqual(data_hollow, True)
        del data_hollow
        self.dm_3x3._validate_shape(data_good)
        self.dm_3x3_cond._validate_shape(data_good)
        del data_good

        # _validate_shap checks just the shape, not the content
        bad_data = np.array([[1.0, 2.0], [3.0, 4.0]])
        data_sym, data_hollow = is_symmetric_and_hollow(bad_data)
        self.assertEqual(data_sym, False)
        del data_sym
        self.assertEqual(data_hollow, False)
        del data_hollow
        data_sym = skbio.stats.distance._utils.is_symmetric(bad_data)
        self.assertEqual(data_sym, False)
        del data_sym
        data_hollow = skbio.stats.distance._utils.is_hollow(bad_data)
        self.assertEqual(data_hollow, False)
        del data_hollow
        self.dm_3x3._validate_shape(bad_data)
        self.dm_3x3_cond._validate_shape(bad_data)
        del bad_data

        # re-try with partially bad data
        bad_data = np.array([[0.0, 2.0], [3.0, 0.0]])
        data_sym, data_hollow = is_symmetric_and_hollow(bad_data)
        self.assertEqual(data_sym, False)
        del data_sym
        self.assertEqual(data_hollow, True)
        del data_hollow
        data_sym = skbio.stats.distance._utils.is_symmetric(bad_data)
        self.assertEqual(data_sym, False)
        del data_sym
        data_hollow = skbio.stats.distance._utils.is_hollow(bad_data)
        self.assertEqual(data_hollow, True)
        del data_hollow
        self.dm_3x3._validate_shape(bad_data)
        self.dm_3x3_cond._validate_shape(bad_data)
        del bad_data

    def test_rename(self):
        # Test successful renaming with a dictionary in strict mode (default)
        dm = DistanceMatrix([[0, 1], [1, 0]], ids=["a", "b"])
        rename_dict = {"a": "x", "b": "y"}
        dm.rename(rename_dict)
        exp = ("x", "y")
        self.assertEqual(dm.ids, exp)

        # Test successful renaming with a function in strict mode (default)
        dm = DistanceMatrix([[0, 1], [1, 0]], ids=["a", "b"])

        def rename_func(x):
            return x + "_1"

        dm.rename(rename_func)
        exp = ("a_1", "b_1")
        self.assertEqual(dm.ids, exp)

        # Test renaming in non-strict mode where one ID is not included in
        # the dictionary
        dm = DistanceMatrix([[0, 1], [1, 0]], ids=["a", "b"])
        rename_dict = {"a": "x"}  # 'b' will retain its original ID
        dm.rename(rename_dict, strict=False)
        exp = ("x", "b")
        self.assertEqual(dm.ids, exp)

        # Test that renaming with strict=True raises an error if not all IDs
        # are included
        dm = DistanceMatrix([[0, 1], [1, 0]], ids=["a", "b"])
        rename_dict = {"a": "x"}  # Missing 'b'
        with self.assertRaises(ValueError):
            dm.rename(rename_dict, strict=True)

    def test_rename_condensed(self):
        # Test successful renaming with a dictionary in strict mode (default)
        dm = DistanceMatrix([[0, 1], [1, 0]], ids=["a", "b"], condensed=True)
        rename_dict = {"a": "x", "b": "y"}
        dm.rename(rename_dict)
        exp = ("x", "y")
        self.assertEqual(dm.ids, exp)

        # Test successful renaming with a function in strict mode (default)
        dm = DistanceMatrix([[0, 1], [1, 0]], ids=["a", "b"], condensed=True)

        def rename_func(x):
            return x + "_1"

        dm.rename(rename_func)
        exp = ("a_1", "b_1")
        self.assertEqual(dm.ids, exp)

        # Test renaming in non-strict mode where one ID is not included in
        # the dictionary
        dm = DistanceMatrix([[0, 1], [1, 0]], ids=["a", "b"], condensed=True)
        rename_dict = {"a": "x"}  # 'b' will retain its original ID
        dm.rename(rename_dict, strict=False)
        exp = ("x", "b")
        self.assertEqual(dm.ids, exp)

        # Test that renaming with strict=True raises an error if not all IDs
        # are included
        dm = DistanceMatrix([[0, 1], [1, 0]], ids=["a", "b"], condensed=True)
        rename_dict = {"a": "x"}  # Missing 'b'
        with self.assertRaises(ValueError):
            dm.rename(rename_dict, strict=True)


class RandomDistanceMatrixTests(TestCase):
    def test_default_usage(self):
        exp = DistanceMatrix(np.asarray([[0.0]]), ["1"])
        obs = randdm(1)
        self.assertEqual(obs, exp)

        obs = randdm(2)
        self.assertEqual(obs.shape, (2, 2))
        self.assertEqual(obs.ids, ("1", "2"))

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
        ids = ["foo", "bar", "baz"]
        obs = randdm(3, ids=ids)
        self.assertEqual(obs.shape, (3, 3))
        self.assertEqual(obs.ids, tuple(ids))

    def test_constructor(self):
        exp = PairwiseMatrix(np.asarray([[0.0]]), ["1"])
        obs = randdm(1, constructor=PairwiseMatrix)
        self.assertEqual(obs, exp)
        self.assertEqual(type(obs), PairwiseMatrix)

    def test_random_fn(self):
        def myrand(size):
            # One dm to rule them all...
            data = np.empty(size)
            data.fill(42)
            return data

        exp = DistanceMatrix(
            np.asarray([[0, 42, 42], [42, 0, 42], [42, 42, 0]]), ["1", "2", "3"]
        )
        obs = randdm(3, random_fn=myrand)
        self.assertEqual(obs, exp)

    def test_random_seed(self):
        obs = randdm(5, random_fn=42).data
        exp = np.array(
            [
                [0.0, 0.97562235, 0.37079802, 0.22723872, 0.75808774],
                [0.97562235, 0.0, 0.92676499, 0.55458479, 0.35452597],
                [0.37079802, 0.92676499, 0.0, 0.06381726, 0.97069802],
                [0.22723872, 0.55458479, 0.06381726, 0.0, 0.89312112],
                [0.75808774, 0.35452597, 0.97069802, 0.89312112, 0.0],
            ]
        )
        npt.assert_almost_equal(obs, exp)

    def test_invalid_input(self):
        # Invalid dimensions.
        with self.assertRaises(PairwiseMatrixError):
            randdm(0)

        # Invalid dimensions.
        with self.assertRaises(ValueError):
            randdm(-1)

        # Invalid number of IDs.
        with self.assertRaises(PairwiseMatrixError):
            randdm(2, ids=["foo"])


class CategoricalStatsHelperFunctionTests(TestCase):
    def setUp(self):
        self.dm = DistanceMatrix(
            [[0.0, 1.0, 2.0], [1.0, 0.0, 3.0], [2.0, 3.0, 0.0]], ["a", "b", "c"]
        )
        self.grouping = [1, 2, 1]
        # Ordering of IDs shouldn't matter, nor should extra IDs.
        self.df = pd.read_csv(
            io.StringIO("ID,Group\nb,Group2\na,Group1\nc,Group1\nd,Group3"), index_col=0
        )
        self.df_missing_id = pd.read_csv(
            io.StringIO("ID,Group\nb,Group2\nc,Group1"), index_col=0
        )

    def test_preprocess_input_with_valid_input(self):
        # Should obtain same result using grouping vector or data frame.
        exp = (
            3,
            2,
            np.array([0, 1, 0]),
            (np.array([0, 0, 1]), np.array([1, 2, 2])),
            np.array([1.0, 2.0, 3.0]),
        )

        obs = _preprocess_input(self.dm, self.grouping, None)
        npt.assert_equal(obs, exp)

        obs = _preprocess_input(self.dm, self.df, "Group")
        npt.assert_equal(obs, exp)

    def test_preprocess_input_raises_error(self):
        # Requires a DistanceMatrix.
        with self.assertRaises(TypeError):
            _preprocess_input(
                PairwiseMatrix([[0, 2], [3, 0]], ["a", "b"]), [1, 2], None
            )

        # Requires column if DataFrame.
        with self.assertRaises(ValueError):
            _preprocess_input(self.dm, self.df, None)

        # Cannot provide column if not data frame.
        with self.assertRaises(ValueError):
            _preprocess_input(self.dm, self.grouping, "Group")

        # Column must exist in data frame.
        with self.assertRaises(ValueError):
            _preprocess_input(self.dm, self.df, "foo")

        # All distance matrix IDs must be in data frame.
        with self.assertRaises(ValueError):
            _preprocess_input(self.dm, self.df_missing_id, "Group")

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


class PairwiseMatrixTests(PairwiseMatrixTestBase, TestCase):
    @classmethod
    def get_matrix_class(cls):
        return PairwiseMatrix

    def setUp(self):
        super(PairwiseMatrixTests, self).setUp()


class SymmetricMatrixTests(SymmetricMatrixTestBase, TestCase):
    @classmethod
    def get_matrix_class(cls):
        return SymmetricMatrix

    def setUp(self):
        super(SymmetricMatrixTests, self).setUp()


class DistanceMatrixTests(DistanceMatrixTestBase, TestCase):
    @classmethod
    def get_matrix_class(cls):
        return DistanceMatrix

    def setUp(self):
        super(DistanceMatrixTests, self).setUp()


class SparseDistanceMatrixTests(TestCase):
    """Test suite for sparse DistanceMatrix implementation."""

    def setUp(self):
        """Set up test data for sparse distance matrices."""
        # Dense test data
        self.dm_3x3_data = [[0.0, 0.01, 4.2], [0.01, 0.0, 12.0], [4.2, 12.0, 0.0]]
        self.dm_3x3_ids = ["a", "b", "c"]

        # Sparse test data (mostly zeros)
        self.sparse_4x4_data = [
            [0.0, 1.0, 0.0, 0.0],
            [1.0, 0.0, 2.0, 0.0],
            [0.0, 2.0, 0.0, 3.0],
            [0.0, 0.0, 3.0, 0.0],
        ]
        self.sparse_4x4_ids = ["w", "x", "y", "z"]

        # Create scipy sparse array
        import scipy.sparse
        self.csr_data = scipy.sparse.csr_array(self.sparse_4x4_data)

    def test_construction_with_sparse_true(self):
        """Test creating sparse distance matrix with sparse=True."""
        dm = DistanceMatrix(self.dm_3x3_data, self.dm_3x3_ids, sparse=True)

        self.assertTrue(dm._flags["SPARSE"])
        self.assertFalse(dm._flags["CONDENSED"])
        self.assertIsInstance(dm.data, scipy.sparse.csr_array)
        self.assertEqual(dm.shape, (3, 3))
        npt.assert_array_almost_equal(dm.data.toarray(), self.dm_3x3_data)

    def test_construction_from_sparse_array(self):
        """Test creating distance matrix from scipy sparse array."""
        dm = DistanceMatrix(self.csr_data, self.sparse_4x4_ids, sparse=True)

        self.assertTrue(dm._flags["SPARSE"])
        self.assertIsInstance(dm.data, scipy.sparse.csr_array)
        npt.assert_array_almost_equal(dm.data.toarray(), self.sparse_4x4_data)

    def test_construction_sparse_false_converts_to_dense(self):
        """Test that sparse=False converts sparse input to dense."""
        dm = DistanceMatrix(self.csr_data, self.sparse_4x4_ids, sparse=False)

        self.assertFalse(dm._flags["SPARSE"])
        self.assertIsInstance(dm.data, np.ndarray)
        npt.assert_array_almost_equal(dm.data, self.sparse_4x4_data)

    def test_construction_sparse_and_condensed_raises_error(self):
        """Test that sparse=True and condensed=True raises error."""
        with self.assertRaisesRegex(
            DistanceMatrixError, "Cannot create sparse matrix in condensed form"
        ):
            DistanceMatrix(
                self.dm_3x3_data, self.dm_3x3_ids, sparse=True, condensed=True
            )

    def test_validation_sparse_symmetry(self):
        """Test that asymmetric sparse matrices are rejected."""
        asym_data = [[0.0, 1.0], [2.0, 0.0]]  # Not symmetric

        with self.assertRaisesRegex(
            DistanceMatrixError,
            "Data must be symmetric"
        ):
            DistanceMatrix(asym_data, ["a", "b"], sparse=True)

    def test_validation_sparse_hollowness(self):
        """Test that non-hollow sparse matrices are rejected."""
        non_hollow = [[1.0, 0.0], [0.0, 1.0]]  # Non-zero diagonal

        with self.assertRaisesRegex(
            DistanceMatrixError,
            "Data must.*be hollow"
        ):
            DistanceMatrix(non_hollow, ["a", "b"], sparse=True)

    def test_validation_sparse_shape(self):
        """Test that non-square sparse matrices are rejected."""
        import scipy.sparse
        non_square = scipy.sparse.csr_array([[0.0, 1.0, 2.0], [1.0, 0.0, 3.0]])

        with self.assertRaisesRegex(
            SymmetricMatrixError,
            "Data must be square"
        ):
            DistanceMatrix(non_square, ["a", "b"], sparse=True)

    def test_shape_property(self):
        """Test that shape property works for sparse matrices."""
        dm = DistanceMatrix(self.sparse_4x4_data, self.sparse_4x4_ids, sparse=True)
        self.assertEqual(dm.shape, (4, 4))

    def test_dtype_property(self):
        """Test that dtype property works for sparse matrices."""
        dm = DistanceMatrix(self.dm_3x3_data, self.dm_3x3_ids, sparse=True)
        self.assertEqual(dm.dtype, np.float64)

    def test_ids_property(self):
        """Test that ids property works for sparse matrices."""
        dm = DistanceMatrix(self.sparse_4x4_data, self.sparse_4x4_ids, sparse=True)
        self.assertEqual(dm.ids, tuple(self.sparse_4x4_ids))

    def test_getitem_single_id(self):
        """Test indexing by single ID returns dense array."""
        dm = DistanceMatrix(self.sparse_4x4_data, self.sparse_4x4_ids, sparse=True)

        row = dm["x"]
        self.assertIsInstance(row, np.ndarray)
        npt.assert_array_almost_equal(row, [1.0, 0.0, 2.0, 0.0])

    def test_getitem_id_pair(self):
        """Test indexing by ID pair returns scalar."""
        dm = DistanceMatrix(self.sparse_4x4_data, self.sparse_4x4_ids, sparse=True)

        val = dm["x", "y"]
        self.assertEqual(val, 2.0)

        val = dm["w", "z"]
        self.assertEqual(val, 0.0)

    def test_getitem_numpy_slicing(self):
        """Test NumPy-style slicing on sparse matrices."""
        dm = DistanceMatrix(self.sparse_4x4_data, self.sparse_4x4_ids, sparse=True)

        # Slicing returns sparse submatrix
        submatrix = dm[1:3, 1:3]
        self.assertIsInstance(submatrix, scipy.sparse.csr_array)
        expected = [[0.0, 2.0], [2.0, 0.0]]
        npt.assert_array_almost_equal(submatrix.toarray(), expected)

    def test_filter_preserves_sparsity(self):
        """Test that filter() preserves sparse format."""
        dm = DistanceMatrix(self.sparse_4x4_data, self.sparse_4x4_ids, sparse=True)

        filtered = dm.filter(["x", "z"])

        self.assertTrue(filtered._flags["SPARSE"])
        self.assertIsInstance(filtered.data, scipy.sparse.csr_array)
        self.assertEqual(filtered.shape, (2, 2))
        self.assertEqual(filtered.ids, ("x", "z"))
        npt.assert_array_almost_equal(
            filtered.data.toarray(),
            [[0.0, 0.0], [0.0, 0.0]]
        )

    def test_copy_preserves_sparsity(self):
        """Test that copy() preserves sparse format."""
        dm = DistanceMatrix(self.sparse_4x4_data, self.sparse_4x4_ids, sparse=True)

        dm_copy = dm.copy()

        self.assertTrue(dm_copy._flags["SPARSE"])
        self.assertIsInstance(dm_copy.data, scipy.sparse.csr_array)
        self.assertIsNot(dm_copy.data, dm.data)
        npt.assert_array_almost_equal(dm_copy.data.toarray(), dm.data.toarray())

    def test_copy_to_condensed_converts_to_dense(self):
        """Test that copying to condensed form densifies sparse matrix."""
        dm = DistanceMatrix(self.sparse_4x4_data, self.sparse_4x4_ids, sparse=True)

        dm_condensed = dm._copy(condensed=True)

        self.assertTrue(dm_condensed._flags["CONDENSED"])
        self.assertFalse(dm_condensed._flags["SPARSE"])
        self.assertIsInstance(dm_condensed.data, np.ndarray)
        self.assertEqual(dm_condensed.data.ndim, 1)

    def test_permute_preserves_sparsity(self):
        """Test that permute() preserves sparse format."""
        dm = DistanceMatrix(self.sparse_4x4_data, self.sparse_4x4_ids, sparse=True)

        dm_permuted = dm.permute(seed=42)

        self.assertTrue(dm_permuted._flags["SPARSE"])
        self.assertIsInstance(dm_permuted.data, scipy.sparse.csr_array)
        self.assertEqual(dm_permuted.shape, (4, 4))

    def test_permute_to_condensed_densifies(self):
        """Test that permute(condensed=True) converts to dense."""
        dm = DistanceMatrix(self.sparse_4x4_data, self.sparse_4x4_ids, sparse=True)

        condensed = dm.permute(condensed=True, seed=42)

        self.assertIsInstance(condensed, np.ndarray)
        self.assertEqual(condensed.ndim, 1)

    def test_redundant_form_converts_to_dense(self):
        """Test that redundant_form() returns dense array."""
        dm = DistanceMatrix(self.sparse_4x4_data, self.sparse_4x4_ids, sparse=True)

        dense = dm.redundant_form()

        self.assertIsInstance(dense, np.ndarray)
        npt.assert_array_almost_equal(dense, self.sparse_4x4_data)

    def test_condensed_form_converts_to_dense(self):
        """Test that condensed_form() returns dense 1D array."""
        dm = DistanceMatrix(self.dm_3x3_data, self.dm_3x3_ids, sparse=True)

        condensed = dm.condensed_form()

        self.assertIsInstance(condensed, np.ndarray)
        self.assertEqual(condensed.ndim, 1)
        # Verify it matches scipy's condensed form
        expected = scipy.spatial.distance.squareform(
            self.dm_3x3_data, force="tovector", checks=False
        )
        npt.assert_array_almost_equal(condensed, expected)

    def test_as_sparse_from_dense(self):
        """Test converting dense matrix to sparse."""
        dm_dense = DistanceMatrix(
            self.sparse_4x4_data, self.sparse_4x4_ids, sparse=False
        )

        dm_sparse = dm_dense.as_sparse()

        self.assertTrue(dm_sparse._flags["SPARSE"])
        self.assertIsInstance(dm_sparse.data, scipy.sparse.csr_array)
        npt.assert_array_almost_equal(
            dm_sparse.data.toarray(),
            self.sparse_4x4_data
        )

    def test_as_sparse_from_sparse_returns_copy(self):
        """Test that as_sparse() on sparse matrix returns copy."""
        dm = DistanceMatrix(self.sparse_4x4_data, self.sparse_4x4_ids, sparse=True)

        dm_sparse = dm.as_sparse()

        self.assertTrue(dm_sparse._flags["SPARSE"])
        self.assertIsNot(dm_sparse, dm)
        self.assertIsNot(dm_sparse.data, dm.data)

    def test_as_dense_from_sparse(self):
        """Test converting sparse matrix to dense."""
        dm_sparse = DistanceMatrix(
            self.sparse_4x4_data, self.sparse_4x4_ids, sparse=True
        )

        dm_dense = dm_sparse.as_dense()

        self.assertFalse(dm_dense._flags["SPARSE"])
        self.assertIsInstance(dm_dense.data, np.ndarray)
        npt.assert_array_almost_equal(dm_dense.data, self.sparse_4x4_data)

    def test_as_dense_from_dense_returns_copy(self):
        """Test that as_dense() on dense matrix returns copy."""
        dm = DistanceMatrix(self.sparse_4x4_data, self.sparse_4x4_ids, sparse=False)

        dm_dense = dm.as_dense()

        self.assertFalse(dm_dense._flags["SPARSE"])
        self.assertIsNot(dm_dense, dm)
        self.assertIsNot(dm_dense.data, dm.data)

    def test_to_data_frame_from_sparse(self):
        """Test conversion to DataFrame from sparse matrix."""
        dm = DistanceMatrix(self.sparse_4x4_data, self.sparse_4x4_ids, sparse=True)

        df = dm.to_data_frame()

        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(df.shape, (4, 4))
        npt.assert_array_almost_equal(df.values, self.sparse_4x4_data)

    def test_sparse_dense_roundtrip(self):
        """Test roundtrip conversion between sparse and dense."""
        dm_original = DistanceMatrix(
            self.sparse_4x4_data, self.sparse_4x4_ids, sparse=False
        )

        dm_sparse = dm_original.as_sparse()
        dm_back = dm_sparse.as_dense()

        npt.assert_array_almost_equal(dm_back.data, dm_original.data)
        self.assertEqual(dm_back.ids, dm_original.ids)


class SparseDistanceMatrixIOTests(TestCase):
    """Test I/O operations for sparse distance matrices."""

    def setUp(self):
        """Set up test data for I/O tests."""
        self.sparse_data = [
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 2.0],
            [0.0, 2.0, 0.0],
        ]
        self.ids = ["a", "b", "c"]

    def test_write_lsmat_densifies_sparse(self):
        """Test that writing sparse matrix to lsmat format densifies it."""
        dm = DistanceMatrix(self.sparse_data, self.ids, sparse=True)

        fh = io.StringIO()
        dm.write(fh, format="lsmat")
        fh.seek(0)

        # Read back as dense
        dm_read = DistanceMatrix.read(fh, format="lsmat")

        self.assertFalse(dm_read._flags["SPARSE"])
        npt.assert_array_almost_equal(dm_read.data, self.sparse_data)

    def test_read_lsmat_with_sparse_true(self):
        """Test reading lsmat format with sparse=True."""
        # Write dense matrix
        dm_dense = DistanceMatrix(self.sparse_data, self.ids, sparse=False)
        fh = io.StringIO()
        dm_dense.write(fh, format="lsmat")
        fh.seek(0)

        # Read as sparse
        dm_sparse = DistanceMatrix.read(fh, format="lsmat", sparse=True)

        self.assertTrue(dm_sparse._flags["SPARSE"])
        self.assertIsInstance(dm_sparse.data, scipy.sparse.csr_array)
        npt.assert_array_almost_equal(dm_sparse.data.toarray(), self.sparse_data)

    def test_lsmat_roundtrip_preserves_values(self):
        """Test that sparse matrix values survive I/O roundtrip."""
        dm_original = DistanceMatrix(self.sparse_data, self.ids, sparse=True)

        # Write and read back
        fh = io.StringIO()
        dm_original.write(fh, format="lsmat")
        fh.seek(0)
        dm_read = DistanceMatrix.read(fh, format="lsmat", sparse=True)

        npt.assert_array_almost_equal(
            dm_read.data.toarray(),
            dm_original.data.toarray()
        )
        self.assertEqual(dm_read.ids, dm_original.ids)


if __name__ == "__main__":
    main()
