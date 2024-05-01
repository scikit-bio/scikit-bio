# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import warnings

import numpy as np
import numpy.testing as npt

from skbio.stats import subsample_counts, isubsample


def setup():
    """Ignore warnings during tests."""
    warnings.simplefilter("ignore")


def teardown():
    """Clear the list of warning filters, so that no filters are active."""
    warnings.resetwarnings()


class SubsampleCountsTests(unittest.TestCase):
    def test_subsample_counts_nonrandom(self):
        a = np.array([0, 5, 0])

        # Subsample same number of items that are in input (without
        # replacement).
        npt.assert_equal(subsample_counts(a, 5), a)

        # Can only choose from one bin.
        exp = np.array([0, 2, 0])
        npt.assert_equal(subsample_counts(a, 2), exp)
        npt.assert_equal(
            subsample_counts(a, 2, replace=True), exp)

        # Subsample zero items.
        a = [3, 0, 1]
        exp = np.array([0, 0, 0])
        npt.assert_equal(subsample_counts(a, 0), exp)
        npt.assert_equal(subsample_counts(a, 0, replace=True), exp)

    def test_subsample_counts_without_replacement(self):
        # Selecting 2 counts from the vector 1000 times yields each of the two
        # possible results at least once each.
        a = np.array([2, 0, 1])
        actual = set()
        for i in range(1000):
            obs = subsample_counts(a, 2)
            actual.add(tuple(obs))
        self.assertEqual(actual, {(1, 0, 1), (2, 0, 0)})

        obs = subsample_counts(a, 2)
        self.assertTrue(np.array_equal(obs, np.array([1, 0, 1])) or
                        np.array_equal(obs, np.array([2, 0, 0])))

    def test_subsample_counts_with_replacement(self):
        # Can choose from all in first bin, all in last bin (since we're
        # sampling with replacement), or split across bins.
        a = np.array([2, 0, 1])
        actual = set()
        for i in range(1000):
            obs = subsample_counts(a, 2, replace=True)
            actual.add(tuple(obs))
        self.assertEqual(actual, {(1, 0, 1), (2, 0, 0), (0, 0, 2)})

        # Test that selecting 35 counts from a 36-count vector 1000 times
        # yields more than 10 different subsamples. If we were subsampling
        # *without* replacement, there would be only 10 possible subsamples
        # because there are 10 nonzero bins in array a. However, there are more
        # than 10 possibilities when sampling *with* replacement.
        a = np.array([2, 0, 1, 2, 1, 8, 6, 0, 3, 3, 5, 0, 0, 0, 5])
        actual = set()
        for i in range(1000):
            obs = subsample_counts(a, 35, replace=True)
            self.assertEqual(obs.sum(), 35)
            actual.add(tuple(obs))
        self.assertTrue(len(actual) > 10)

    def test_subsample_counts_with_replacement_equal_n(self):
        # test when n == counts.sum()
        a = np.array([0, 0, 3, 4, 2, 1])
        actual = set()
        for i in range(1000):
            obs = subsample_counts(a, 10, replace=True)
            self.assertEqual(obs.sum(), 10)
            actual.add(tuple(obs))
        self.assertTrue(len(actual) > 1)

    def test_subsample_counts_invalid_input(self):
        # Negative n.
        with self.assertRaises(ValueError):
            subsample_counts([1, 2, 3], -1)

        # Wrong number of dimensions.
        with self.assertRaises(ValueError):
            subsample_counts([[1, 2, 3], [4, 5, 6]], 2)

        # Input has too few counts.
        with self.assertRaises(ValueError):
            subsample_counts([0, 5, 0], 6, replace=False)

        # Input has too counts, but should work with bootstrap
        subsample_counts([0, 5, 0], 6, replace=True)


class ISubsampleTests(unittest.TestCase):
    def setUp(self):
        np.random.seed(123)

        # comment indicates the expected random value
        self.sequences = [
            ('a_1', 'AATTGGCC-a1'),  # 2, 3624216819017203053
            ('a_2', 'AATTGGCC-a2'),  # 5, 5278339153051796802
            ('b_1', 'AATTGGCC-b1'),  # 4, 4184670734919783522
            ('b_2', 'AATTGGCC-b2'),  # 0, 946590342492863505
            ('a_4', 'AATTGGCC-a4'),  # 3, 4048487933969823850
            ('a_3', 'AATTGGCC-a3'),  # 7, 7804936597957240377
            ('c_1', 'AATTGGCC-c1'),  # 8, 8868534167180302049
            ('a_5', 'AATTGGCC-a5'),  # 1, 3409506807702804593
            ('c_2', 'AATTGGCC-c2'),  # 9, 8871627813779918895
            ('c_3', 'AATTGGCC-c3')   # 6, 7233291490207274528
        ]

    def mock_sequence_iter(self, items):
        return ({'SequenceID': sid, 'Sequence': seq} for sid, seq in items)

    def test_isubsample_simple(self):
        maximum = 10

        def bin_f(x):
            return x['SequenceID'].rsplit('_', 1)[0]

        # note, the result here is sorted by sequence_id but is in heap order
        # by the random values associated to each sequence
        exp = sorted([('a', {'SequenceID': 'a_5', 'Sequence': 'AATTGGCC-a5'}),
                      ('a', {'SequenceID': 'a_1', 'Sequence': 'AATTGGCC-a1'}),
                      ('a', {'SequenceID': 'a_4', 'Sequence': 'AATTGGCC-a4'}),
                      ('a', {'SequenceID': 'a_3', 'Sequence': 'AATTGGCC-a3'}),
                      ('a', {'SequenceID': 'a_2', 'Sequence': 'AATTGGCC-a2'}),
                      ('b', {'SequenceID': 'b_2', 'Sequence': 'AATTGGCC-b2'}),
                      ('b', {'SequenceID': 'b_1', 'Sequence': 'AATTGGCC-b1'}),
                      ('c', {'SequenceID': 'c_3', 'Sequence': 'AATTGGCC-c3'}),
                      ('c', {'SequenceID': 'c_2', 'Sequence': 'AATTGGCC-c2'}),
                      ('c', {'SequenceID': 'c_1', 'Sequence': 'AATTGGCC-c1'})],
                     key=lambda x: x[0])
        obs = isubsample(self.mock_sequence_iter(self.sequences), maximum,
                         bin_f=bin_f)
        self.assertEqual(sorted(obs, key=lambda x: x[0]), exp)

    def test_per_sample_sequences_min_seqs(self):
        maximum = 10
        minimum = 3

        def bin_f(x):
            return x['SequenceID'].rsplit('_', 1)[0]

        # note, the result here is sorted by sequence_id but is in heap order
        # by the random values associated to each sequence
        exp = sorted([('a', {'SequenceID': 'a_5', 'Sequence': 'AATTGGCC-a5'}),
                      ('a', {'SequenceID': 'a_1', 'Sequence': 'AATTGGCC-a1'}),
                      ('a', {'SequenceID': 'a_4', 'Sequence': 'AATTGGCC-a4'}),
                      ('a', {'SequenceID': 'a_3', 'Sequence': 'AATTGGCC-a3'}),
                      ('a', {'SequenceID': 'a_2', 'Sequence': 'AATTGGCC-a2'}),
                      ('c', {'SequenceID': 'c_3', 'Sequence': 'AATTGGCC-c3'}),
                      ('c', {'SequenceID': 'c_2', 'Sequence': 'AATTGGCC-c2'}),
                      ('c', {'SequenceID': 'c_1', 'Sequence': 'AATTGGCC-c1'})],
                     key=lambda x: x[0])
        obs = isubsample(self.mock_sequence_iter(self.sequences), maximum,
                         minimum, bin_f=bin_f)
        self.assertEqual(sorted(obs, key=lambda x: x[0]), exp)

    def test_per_sample_sequences_complex(self):
        maximum = 2

        def bin_f(x):
            return x['SequenceID'].rsplit('_', 1)[0]

        exp = sorted([('a', {'SequenceID': 'a_2', 'Sequence': 'AATTGGCC-a2'}),
                      ('a', {'SequenceID': 'a_3', 'Sequence': 'AATTGGCC-a3'}),
                      ('b', {'SequenceID': 'b_2', 'Sequence': 'AATTGGCC-b2'}),
                      ('b', {'SequenceID': 'b_1', 'Sequence': 'AATTGGCC-b1'}),
                      ('c', {'SequenceID': 'c_1', 'Sequence': 'AATTGGCC-c1'}),
                      ('c', {'SequenceID': 'c_2', 'Sequence': 'AATTGGCC-c2'})],
                     key=lambda x: x[0])
        obs = isubsample(self.mock_sequence_iter(self.sequences), maximum,
                         bin_f=bin_f, buf_size=1)
        self.assertEqual(sorted(obs, key=lambda x: x[0]), exp)

    def test_min_gt_max(self):
        gen = isubsample([1, 2, 3], maximum=2, minimum=10)
        with self.assertRaises(ValueError):
            next(gen)

    def test_min_lt_zero(self):
        gen = isubsample([1, 2, 3], maximum=0, minimum=-10)
        with self.assertRaises(ValueError):
            next(gen)

    def test_max_lt_zero(self):
        gen = isubsample([1, 2, 3], maximum=-10)
        with self.assertRaises(ValueError):
            next(gen)

    def test_binf_is_none(self):
        maximum = 2
        items = [1, 2]
        exp = [(True, 1), (True, 2)]
        obs = isubsample(items, maximum)
        self.assertEqual(list(obs), exp)


if __name__ == '__main__':
    unittest.main()
