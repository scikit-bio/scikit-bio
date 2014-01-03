#!/usr/bin/env python
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from itertools import izip
from StringIO import StringIO

import numpy as np

from bipy.core.distance import (DistanceMatrix, DistanceMatrixError,
        DistanceMatrixFormatError, MissingDataError, MissingHeaderError,
        MissingSampleIDError, SampleIDMismatchError)
from bipy.util.unit_test import TestCase, main

class DistanceMatrixTests(TestCase):
    """Tests for the DistanceMatrix class."""

    def setUp(self):
        """Set up test data for use in distance matrix unit tests."""
        dm_1x1_data = [[0.0]]
        self.dm_1x1_f = StringIO(DM_1x1_F)
        self.dm_1x1 = DistanceMatrix(dm_1x1_data, ['a'])

        dm_2x2_data = [[0.0, 0.123], [0.123, 0.0]]
        self.dm_2x2_f = StringIO(DM_2x2_F)
        self.dm_2x2 = DistanceMatrix(dm_2x2_data, ['a', 'b'])

        dm_3x3_data = [[0.0, 0.01, 4.2], [0.01, 0.0, 12.0], [4.2, 12.0, 0.0]]
        self.dm_3x3_f = StringIO(DM_3x3_F)
        self.dm_3x3 = DistanceMatrix(dm_3x3_data, ['a', 'b', 'c'])

        self.dm_3x3_whitespace_f = StringIO(DM_3x3_WHITESPACE_F)

        self.dm_f_lines = [DM_1x1_F, DM_2x2_F, DM_3x3_F]
        self.dm_fs = [self.dm_1x1_f, self.dm_2x2_f, self.dm_3x3_f]
        self.dms = [self.dm_1x1, self.dm_2x2, self.dm_3x3]
        self.dm_shapes = [(1, 1), (2, 2), (3, 3)]
        self.dm_sizes = [1, 4, 9]
        self.dm_condensed_forms = [np.array([]), np.array([0.123]),
                                   np.array([0.01, 4.2, 12.0])]
        self.dm_redundant_forms = [np.array(dm_1x1_data),
                                   np.array(dm_2x2_data),
                                   np.array(dm_3x3_data)]

        self.bad_dm_f1 = StringIO(BAD_DM_F1)
        self.bad_dm_f2 = StringIO(BAD_DM_F2)
        self.bad_dm_f3 = StringIO(BAD_DM_F3)
        self.bad_dm_f4 = StringIO(BAD_DM_F4)
        self.bad_dm_f5 = StringIO(BAD_DM_F5)
        self.bad_dm_f6 = StringIO(BAD_DM_F6)
        self.bad_dm_f7 = StringIO(BAD_DM_F7)

    def test_round_trip_read_write(self):
        """Test reading, writing, and reading again works as expected."""
        for dm_f in self.dm_fs:
            # Read.
            dm1 = DistanceMatrix.from_file(dm_f)

            # Write.
            out_f = StringIO()
            dm1.to_file(out_f)
            out_f.seek(0)

            # Read.
            dm2 = DistanceMatrix.from_file(out_f)
            self.assertEqual(dm1, dm2)

    def test_from_file(self):
        """Should parse and return a valid DistanceMatrix given a file."""
        for dm_f, dm in izip(self.dm_fs, self.dms):
            obs = DistanceMatrix.from_file(dm_f)
            self.assertEqual(obs, dm)

        # Correctly parses file with extra empty (whitespace-only) lines at
        # end.
        obs = DistanceMatrix.from_file(self.dm_3x3_whitespace_f)
        self.assertEqual(obs, self.dm_3x3)

    def test_from_file_invalid_input(self):
        """Raises error on ill-formatted distance matrix file."""
        # Empty dm.
        with self.assertRaises(MissingHeaderError):
            _ = DistanceMatrix.from_file([])

        # Number of values don't match number of sample IDs.
        with self.assertRaises(DistanceMatrixFormatError):
            _ = DistanceMatrix.from_file(self.bad_dm_f1)

        # Mismatched sample IDs.
        with self.assertRaises(SampleIDMismatchError):
            _ = DistanceMatrix.from_file(self.bad_dm_f2)

        # Extra data at end.
        with self.assertRaises(DistanceMatrixFormatError):
            _ = DistanceMatrix.from_file(self.bad_dm_f3)

        # Missing data.
        with self.assertRaises(MissingDataError):
            _ = DistanceMatrix.from_file(self.bad_dm_f4)

        # Header, but no data.
        with self.assertRaises(MissingDataError):
            _ = DistanceMatrix.from_file(self.bad_dm_f5)

        # Nonsymmetric.
        with self.assertRaises(DistanceMatrixError):
            _ = DistanceMatrix.from_file(self.bad_dm_f6)

        # Non-hollow.
        with self.assertRaises(DistanceMatrixError):
            _ = DistanceMatrix.from_file(self.bad_dm_f7)

    def test_to_file(self):
        """Should serialize a DistanceMatrix to file."""
        for dm_f_line, dm in izip(self.dm_f_lines, self.dms):
            for conserve_memory in True, False:
                obs_f = StringIO()
                dm.to_file(obs_f, conserve_memory=conserve_memory)
                obs = obs_f.getvalue()
                obs_f.close()

                self.assertEqual(obs, dm_f_line)

    def test_init_invalid_input(self):
        """Raises error on invalid distance matrix data / sample IDs."""
        # Empty data.
        with self.assertRaises(DistanceMatrixError):
            _ = DistanceMatrix([], [])

        # Invalid number of dimensions.
        with self.assertRaises(DistanceMatrixError):
            _ = DistanceMatrix([1, 2, 3], ['a'])

        # Dimensions don't match.
        with self.assertRaises(DistanceMatrixError):
            _ = DistanceMatrix([[1, 2, 3]], ['a'])

        data = [[0, 1], [1, 0]]

        # Duplicate sample IDs.
        with self.assertRaises(DistanceMatrixError):
            _ = DistanceMatrix(data, ['a', 'a'])

        # Number of sample IDs don't match dimensions.
        with self.assertRaises(DistanceMatrixError):
            _ = DistanceMatrix(data, ['a', 'b', 'c'])

        # Non-hollow.
        data = [[0.0, 1.0], [1.0, 0.01]]
        with self.assertRaises(DistanceMatrixError):
            _ = DistanceMatrix(data, ['a', 'b'])

        # Nonsymmetric.
        data = [[0.0, 2.0], [1.0, 0.0]]
        with self.assertRaises(DistanceMatrixError):
            _ = DistanceMatrix(data, ['a', 'b'])

    def test_data(self):
        """Test retrieving/setting data matrix."""
        for dm, exp in izip(self.dms, self.dm_redundant_forms):
            obs = dm.data
            self.assertTrue(np.array_equal(obs, exp))

        with self.assertRaises(AttributeError):
            self.dm_3x3.data = 'foo'

    def test_sample_ids(self):
        """Test retrieving/setting sample IDs."""
        obs = self.dm_3x3.sample_ids
        self.assertEqual(obs, ('a', 'b', 'c'))

        # Test that we overwrite the existing sample IDs and that the sample
        # index is correctly rebuilt.
        new_sids = ['foo', 'bar', 'baz']
        self.dm_3x3.sample_ids = new_sids
        obs = self.dm_3x3.sample_ids
        self.assertEqual(obs, tuple(new_sids))
        self.assertTrue(np.array_equal(self.dm_3x3['bar'],
                                       np.array([0.01, 0.0, 12.0])))
        with self.assertRaises(MissingSampleIDError):
            _ = self.dm_3x3['b']

    def test_sample_ids_invalid_input(self):
        """Test setting invalid sample IDs raises an error."""
        with self.assertRaises(DistanceMatrixError):
            self.dm_3x3.sample_ids = ['foo', 'bar']
        # Make sure that we can still use the distance matrix after trying to
        # be evil.
        obs = self.dm_3x3.sample_ids
        self.assertEqual(obs, ('a', 'b', 'c'))

    def test_dtype(self):
        """Test retrieving dtype of data matrix."""
        for dm in self.dms:
            self.assertEqual(dm.dtype, np.float64)

    def test_shape(self):
        """Test retrieving shape of data matrix."""
        for dm, shape in izip(self.dms, self.dm_shapes):
            self.assertEqual(dm.shape, shape)

    def test_num_samples(self):
        """Test retrieving the number of samples in the distance matrix."""
        for dm, shape in izip(self.dms, self.dm_shapes):
            self.assertEqual(dm.num_samples, shape[0])

    def test_size(self):
        """Test retrieving size of data matrix."""
        for dm, size in izip(self.dms, self.dm_sizes):
            self.assertEqual(dm.size, size)

    def test_transpose(self):
        """Test retrieving transpose of distance matrix."""
        for dm in self.dms:
            self.assertEqual(dm.T, dm)
            self.assertEqual(dm.transpose(), dm)
            self.assertTrue(dm.transpose() is dm)

    def test_condensed_form(self):
        """Test retrieving the data matrix in condensed form."""
        for dm, condensed in izip(self.dms, self.dm_condensed_forms):
            obs = dm.condensed_form()
            self.assertTrue(np.array_equal(obs, condensed))

    def test_redundant_form(self):
        """Test retrieving the data matrix in redundant form."""
        for dm, redundant in izip(self.dms, self.dm_redundant_forms):
            obs = dm.redundant_form()
            self.assertTrue(np.array_equal(obs, redundant))

    def test_copy(self):
        """Test correct copying of a DistanceMatrix."""
        copy = self.dm_2x2.copy()
        self.assertEqual(copy, self.dm_2x2)
        self.assertFalse(copy.data is self.dm_2x2.data)
        # deepcopy doesn't actually create a copy of the IDs because it is a
        # tuple of strings, which is fully immutable.
        self.assertTrue(copy.sample_ids is self.dm_2x2.sample_ids)

        new_ids = ['hello', 'world']
        copy.sample_ids = new_ids
        self.assertNotEqual(copy.sample_ids, self.dm_2x2.sample_ids)

        copy = self.dm_2x2.copy()
        copy.data[0,1] = 0.0001
        self.assertFalse(np.array_equal(copy.data, self.dm_2x2.data))

    def test_str(self):
        """Test retrieving string representation of a DistanceMatrix."""
        for dm in self.dms:
            obs = str(dm)
            # Do some very light testing here to make sure we're getting a
            # non-empty string back. We don't want to test the exact
            # formatting.
            self.assertTrue(obs)

    def test_eq(self):
        """Test DistanceMatrix equality test functions correctly."""
        for dm in self.dms:
            copy = dm.copy()
            self.assertTrue(dm == dm)
            self.assertTrue(copy == copy)
            self.assertTrue(dm == copy)
            self.assertTrue(copy == dm)

        self.assertFalse(self.dm_1x1 == self.dm_3x3)

    def test_ne(self):
        """Test unequal dms are identified as such."""
        # Wrong class.
        self.assertTrue(self.dm_3x3 != 'foo')

        # Wrong shape.
        self.assertTrue(self.dm_3x3 != self.dm_1x1)

        # Wrong sample IDs.
        other = self.dm_3x3.copy()
        other.sample_ids = ['foo', 'bar', 'baz']
        self.assertTrue(self.dm_3x3 != other)

        # Wrong data.
        other = self.dm_3x3.copy()
        other.data[1,0] = 42.42
        self.assertTrue(self.dm_3x3 != other)

        self.assertFalse(self.dm_2x2 != self.dm_2x2)

    def test_getitem(self):
        """Test retrieving vectors by sample ID."""
        obs = self.dm_1x1['a']
        self.assertTrue(np.array_equal(obs, np.array([0.0])))

        obs = self.dm_3x3['c']
        self.assertTrue(np.array_equal(obs, np.array([4.2, 12.0, 0.0])))

        with self.assertRaises(MissingSampleIDError):
            _ = self.dm_2x2['c']

    def test_validate(self):
        """Empty stub: DistanceMatrix._validate already tested elsewhere."""
        pass

    def test_index_list(self):
        """Empty stub: DistanceMatrix._index_list already tested elsewhere."""
        pass

    def test_format_sample_ids(self):
        """Empty stub: DistanceMatrix._format_sample_ids tested elsewhere."""
        pass

    def test_pprint_sample_ids(self):
        """Test pretty-print formatting of sample IDs."""
        # No truncation.
        exp = 'a, b, c'
        obs = self.dm_3x3._pprint_sample_ids()
        self.assertEqual(obs, exp)

        # Truncation.
        exp = 'a, b, ...'
        obs = self.dm_3x3._pprint_sample_ids(max_chars=5)
        self.assertEqual(obs, exp)

# 1x1:
#     0.0
DM_1x1_F = "\ta\na\t0.0\n"

# 2x2:
#       0.0  0.123
#     0.123    0.0
DM_2x2_F = "\ta\tb\na\t0.0\t0.123\nb\t0.123\t0.0\n"

# 3x3:
#      0.0   0.01   4.2
#     0.01    0.0  12.0
#      4.2   12.0   0.0
DM_3x3_F = ("\ta\tb\tc\na\t0.0\t0.01\t4.2\nb\t0.01\t0.0\t12.0\n"
            "c\t4.2\t12.0\t0.0\n")

# Extra whitespace-only lines at end.
DM_3x3_WHITESPACE_F = ("\ta\tb\tc\na\t0.0\t0.01\t4.2\nb\t0.01\t0.0\t12.0\n"
                       "c\t4.2\t12.0\t0.0\n\n   \t \n\t\t\t\n ")

# missing data
BAD_DM_F1 = 'a\tb\na\t0\t1\nb\t1'

# mismatched sample IDs
BAD_DM_F2 = '\ta\tb\nb\t0\t1\na\t1\t0'

# extra data lines
BAD_DM_F3 = '\ta\tb\na\t0\t1\nb\t1\t0\nfoo'

# missing data lines
BAD_DM_F4 = '\ta\tb\na\t0\t1\n'

# no data lines
BAD_DM_F5 = '\ta\tb\n'

# nonsymmetric
BAD_DM_F6 = '\ta\tb\na\t0\t1\nb\t2\t0\n'

# non-hollow
BAD_DM_F7 = '\ta\tb\na\t0\t1\nb\t1\t0.1\n'

if __name__ == '__main__':
    main()
