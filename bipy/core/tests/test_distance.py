#!/usr/bin/env python
"""Unit tests for the bipy.core.distance module."""
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, bipy development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from itertools import izip
from StringIO import StringIO
from tempfile import TemporaryFile

import numpy as np

from bipy.core.distance import (random_distance_matrix, DistanceMatrix,
                                DistanceMatrixError, DistanceMatrixFormatError,
                                MissingDataError, MissingHeaderError,
                                MissingIDError, IDMismatchError,
                                SymmetricDistanceMatrix)
from bipy.util.unit_test import TestCase, main


class DistanceMatrixTestData(TestCase):
    """Test data used in DistanceMatrix and subclass unit tests."""

    def setUp(self):
        """Set up test data for use in distance matrix unit tests."""
        self.dm_1x1_data = [[0.0]]
        self.dm_1x1_f = StringIO(DM_1x1_F)

        self.dm_2x2_data = [[0.0, 0.123], [0.123, 0.0]]
        self.dm_2x2_f = StringIO(DM_2x2_F)

        self.dm_2x2_asym_data = [[0.0, 1.0], [-2.0, 0.0]]
        self.dm_2x2_asym_f = StringIO(DM_2x2_ASYM_F)

        self.dm_3x3_data = [[0.0, 0.01, 4.2], [0.01, 0.0, 12.0],
                            [4.2, 12.0, 0.0]]
        self.dm_3x3_f = StringIO(DM_3x3_F)

        self.dm_3x3_whitespace_f = StringIO('\n'.join(DM_3x3_WHITESPACE_F))

        self.tmp_f = TemporaryFile(prefix='bipy.core.tests.test_distance',
                                   suffix='.txt')

        self.bad_dm_f1 = StringIO(BAD_DM_F1)
        self.bad_dm_f2 = StringIO(BAD_DM_F2)
        self.bad_dm_f3 = StringIO(BAD_DM_F3)
        self.bad_dm_f4 = StringIO(BAD_DM_F4)
        self.bad_dm_f5 = StringIO(BAD_DM_F5)
        self.bad_dm_f6 = StringIO(BAD_DM_F6)

    def tearDown(self):
        """Delete any temporary files."""
        self.tmp_f.close()


class DistanceMatrixTests(DistanceMatrixTestData):
    """Tests for the DistanceMatrix class."""

    def setUp(self):
        """Create some DistanceMatrix instances for use in tests."""
        super(DistanceMatrixTests, self).setUp()

        self.dm_1x1 = DistanceMatrix(self.dm_1x1_data, ['a'])
        self.dm_2x2 = DistanceMatrix(self.dm_2x2_data, ['a', 'b'])
        self.dm_2x2_asym = DistanceMatrix(self.dm_2x2_asym_data, ['a', 'b'])
        self.dm_3x3 = DistanceMatrix(self.dm_3x3_data, ['a', 'b', 'c'])

        self.dms = [self.dm_1x1, self.dm_2x2, self.dm_2x2_asym, self.dm_3x3]
        self.dm_f_lines = [DM_1x1_F, DM_2x2_F, DM_2x2_ASYM_F, DM_3x3_F]
        self.dm_fs = [self.dm_1x1_f, self.dm_2x2_f, self.dm_2x2_asym_f,
                      self.dm_3x3_f]
        self.dm_shapes = [(1, 1), (2, 2), (2, 2), (3, 3)]
        self.dm_sizes = [1, 4, 4, 9]
        self.dm_transposes = [self.dm_1x1, self.dm_2x2,
                              DistanceMatrix([[0, -2], [1, 0]], ['a', 'b']),
                              self.dm_3x3]
        self.dm_redundant_forms = [np.array(self.dm_1x1_data),
                                   np.array(self.dm_2x2_data),
                                   np.array(self.dm_2x2_asym_data),
                                   np.array(self.dm_3x3_data)]

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

    def test_from_file_extra_junk(self):
        """Should correctly parse a file with extra whitespace and comments."""
        obs = DistanceMatrix.from_file(self.dm_3x3_whitespace_f)
        self.assertEqual(obs, self.dm_3x3)

    def test_from_file_list_of_strings(self):
        """Should correctly parse a list of strings."""
        obs = DistanceMatrix.from_file(DM_3x3_WHITESPACE_F)
        self.assertEqual(obs, self.dm_3x3)

    def test_from_file_real_file(self):
        """Should correctly parse a real on-disk file."""
        self.tmp_f.write('\n'.join(DM_3x3_WHITESPACE_F))
        self.tmp_f.seek(0)

        obs = DistanceMatrix.from_file(self.tmp_f)
        self.assertEqual(obs, self.dm_3x3)

    def test_from_file_invalid_input(self):
        """Raises error on ill-formatted distance matrix file."""
        # Empty dm.
        with self.assertRaises(MissingHeaderError):
            _ = DistanceMatrix.from_file([])

        # Number of values don't match number of IDs.
        with self.assertRaises(DistanceMatrixFormatError):
            _ = DistanceMatrix.from_file(self.bad_dm_f1)

        # Mismatched IDs.
        with self.assertRaises(IDMismatchError):
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

        # Non-hollow.
        with self.assertRaises(DistanceMatrixError):
            _ = DistanceMatrix.from_file(self.bad_dm_f6)

    def test_to_file(self):
        """Should serialize a DistanceMatrix to file."""
        for dm_f_line, dm in izip(self.dm_f_lines, self.dms):
            obs_f = StringIO()
            dm.to_file(obs_f)
            obs = obs_f.getvalue()
            obs_f.close()

            self.assertEqual(obs, dm_f_line)

    def test_init_from_dm(self):
        """Constructs a dm from a dm."""
        ids = ['foo', 'bar', 'baz']

        # DistanceMatrix -> DistanceMatrix
        exp = DistanceMatrix(self.dm_3x3_data, ids)
        obs = DistanceMatrix(self.dm_3x3, ids)
        self.assertEqual(obs, exp)
        # Test that copy of data is not made.
        self.assertTrue(obs.data is self.dm_3x3.data)
        obs.data[0, 1] = 424242
        self.assertTrue(np.array_equal(obs.data, self.dm_3x3.data))

        # SymmetricDistanceMatrix -> DistanceMatrix
        exp = DistanceMatrix(self.dm_3x3_data, ids)
        obs = DistanceMatrix(SymmetricDistanceMatrix(self.dm_3x3_data,
                                                     ('a', 'b', 'c')), ids)
        self.assertEqual(obs, exp)

        # DistanceMatrix -> SymmetricDistanceMatrix
        with self.assertRaises(DistanceMatrixError):
            _ = SymmetricDistanceMatrix(self.dm_2x2_asym, ['foo', 'bar'])

    def test_init_invalid_input(self):
        """Raises error on invalid distance matrix data / IDs."""
        # Empty data.
        with self.assertRaises(DistanceMatrixError):
            _ = DistanceMatrix([], [])

        # Another type of empty data.
        with self.assertRaises(DistanceMatrixError):
            _ = DistanceMatrix(np.empty((0, 0)), [])

        # Invalid number of dimensions.
        with self.assertRaises(DistanceMatrixError):
            _ = DistanceMatrix([1, 2, 3], ['a'])

        # Dimensions don't match.
        with self.assertRaises(DistanceMatrixError):
            _ = DistanceMatrix([[1, 2, 3]], ['a'])

        data = [[0, 1], [1, 0]]

        # Duplicate IDs.
        with self.assertRaises(DistanceMatrixError):
            _ = DistanceMatrix(data, ['a', 'a'])

        # Number of IDs don't match dimensions.
        with self.assertRaises(DistanceMatrixError):
            _ = DistanceMatrix(data, ['a', 'b', 'c'])

        # Non-hollow.
        data = [[0.0, 1.0], [1.0, 0.01]]
        with self.assertRaises(DistanceMatrixError):
            _ = DistanceMatrix(data, ['a', 'b'])

    def test_data(self):
        """Test retrieving/setting data matrix."""
        for dm, exp in izip(self.dms, self.dm_redundant_forms):
            obs = dm.data
            self.assertTrue(np.array_equal(obs, exp))

        with self.assertRaises(AttributeError):
            self.dm_3x3.data = 'foo'

    def test_ids(self):
        """Test retrieving/setting IDs."""
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
            _ = self.dm_3x3['b']

    def test_ids_invalid_input(self):
        """Test setting invalid IDs raises an error."""
        with self.assertRaises(DistanceMatrixError):
            self.dm_3x3.ids = ['foo', 'bar']
        # Make sure that we can still use the distance matrix after trying to
        # be evil.
        obs = self.dm_3x3.ids
        self.assertEqual(obs, ('a', 'b', 'c'))

    def test_dtype(self):
        """Test retrieving dtype of data matrix."""
        for dm in self.dms:
            self.assertEqual(dm.dtype, np.float64)

    def test_shape(self):
        """Test retrieving shape of data matrix."""
        for dm, shape in izip(self.dms, self.dm_shapes):
            self.assertEqual(dm.shape, shape)

    def test_size(self):
        """Test retrieving size of data matrix."""
        for dm, size in izip(self.dms, self.dm_sizes):
            self.assertEqual(dm.size, size)

    def test_transpose(self):
        """Test retrieving transpose of distance matrix."""
        for dm, transpose in izip(self.dms, self.dm_transposes):
            self.assertEqual(dm.T, transpose)
            self.assertEqual(dm.transpose(), transpose)
            # We should get a reference to a different object back, even if the
            # transpose is the same as the original.
            self.assertTrue(dm.transpose() is not dm)

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
        self.assertTrue(copy.ids is self.dm_2x2.ids)

        new_ids = ['hello', 'world']
        copy.ids = new_ids
        self.assertNotEqual(copy.ids, self.dm_2x2.ids)

        copy = self.dm_2x2.copy()
        copy.data[0, 1] = 0.0001
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

        # Wrong IDs.
        other = self.dm_3x3.copy()
        other.ids = ['foo', 'bar', 'baz']
        self.assertTrue(self.dm_3x3 != other)

        # Wrong data.
        other = self.dm_3x3.copy()
        other.data[1, 0] = 42.42
        self.assertTrue(self.dm_3x3 != other)

        self.assertFalse(self.dm_2x2 != self.dm_2x2)

    def test_getslice(self):
        """Test that __getslice__ defers to __getitem__."""
        # Slice of first dimension only.
        obs = self.dm_2x2[1:]
        self.assertTrue(np.array_equal(obs, np.array([[0.123, 0.0]])))
        self.assertEqual(type(obs), np.ndarray)

    def test_getitem_by_id(self):
        """Test retrieving row vectors by ID."""
        obs = self.dm_1x1['a']
        self.assertTrue(np.array_equal(obs, np.array([0.0])))

        obs = self.dm_2x2_asym['b']
        self.assertTrue(np.array_equal(obs, np.array([-2.0, 0.0])))

        obs = self.dm_3x3['c']
        self.assertTrue(np.array_equal(obs, np.array([4.2, 12.0, 0.0])))

        with self.assertRaises(MissingIDError):
            _ = self.dm_2x2['c']

    def test_getitem_by_id_pair(self):
        """Test retrieving elements by ID pair."""
        # Same object.
        self.assertEqual(self.dm_1x1['a', 'a'], 0.0)

        # Different objects (symmetric).
        self.assertEqual(self.dm_3x3['b', 'c'], 12.0)
        self.assertEqual(self.dm_3x3['c', 'b'], 12.0)

        # Different objects (asymmetric).
        self.assertEqual(self.dm_2x2_asym['a', 'b'], 1.0)
        self.assertEqual(self.dm_2x2_asym['b', 'a'], -2.0)

        with self.assertRaises(MissingIDError):
            _ = self.dm_2x2['a', 'c']

    def test_getitem_ndarray_indexing(self):
        """Test __getitem__ delegates to underlying ndarray."""
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
            _ = self.dm_3x3[:, 3]

    def test_parse_ids(self):
        """Empty stub: DistanceMatrix._parse_ids tested elsewhere."""
        pass

    def test_validate(self):
        """Empty stub: DistanceMatrix._validate already tested elsewhere."""
        pass

    def test_index_list(self):
        """Empty stub: DistanceMatrix._index_list already tested elsewhere."""
        pass

    def test_is_id_pair(self):
        """Empty stub: DistanceMatrix._is_id_pair already tested elsewhere."""
        pass

    def test_format_ids(self):
        """Empty stub: DistanceMatrix._format_ids tested elsewhere."""
        pass

    def test_pprint_ids(self):
        """Test pretty-print formatting of IDs."""
        # No truncation.
        exp = 'a, b, c'
        obs = self.dm_3x3._pprint_ids()
        self.assertEqual(obs, exp)

        # Truncation.
        exp = 'a, b, ...'
        obs = self.dm_3x3._pprint_ids(max_chars=5)
        self.assertEqual(obs, exp)


class SymmetricDistanceMatrixTests(DistanceMatrixTestData):
    """Tests for the SymmetricDistanceMatrix class."""

    def setUp(self):
        """Create some SymmetricDistanceMatrix instances for use in tests."""
        super(SymmetricDistanceMatrixTests, self).setUp()

        self.dm_1x1 = SymmetricDistanceMatrix(self.dm_1x1_data, ['a'])
        self.dm_2x2 = SymmetricDistanceMatrix(self.dm_2x2_data, ['a', 'b'])
        self.dm_3x3 = SymmetricDistanceMatrix(self.dm_3x3_data,
                                              ['a', 'b', 'c'])

        self.dms = [self.dm_1x1, self.dm_2x2, self.dm_3x3]
        self.dm_condensed_forms = [np.array([]), np.array([0.123]),
                                   np.array([0.01, 4.2, 12.0])]

    def test_from_file_invalid_input(self):
        """Raises error on invalid distance matrix file."""
        # Asymmetric.
        with self.assertRaises(DistanceMatrixError):
            _ = SymmetricDistanceMatrix.from_file(self.dm_2x2_asym_f)

    def test_init_invalid_input(self):
        """Raises error on invalid distance matrix data / IDs."""
        # Asymmetric.
        data = [[0.0, 2.0], [1.0, 0.0]]
        with self.assertRaises(DistanceMatrixError):
            _ = SymmetricDistanceMatrix(data, ['a', 'b'])

        # Ensure that the superclass validation is still being performed.
        with self.assertRaises(DistanceMatrixError):
            _ = SymmetricDistanceMatrix([[1, 2, 3]], ['a'])

    def test_condensed_form(self):
        """Test retrieving the data matrix in condensed form."""
        for dm, condensed in izip(self.dms, self.dm_condensed_forms):
            obs = dm.condensed_form()
            self.assertTrue(np.array_equal(obs, condensed))

    def test_eq(self):
        """Test data equality between different distance matrix types."""
        # Compare SymmetricDistanceMatrix to DistanceMatrix, where both have
        # the same data and IDs.
        eq_dm = DistanceMatrix(self.dm_3x3_data, ['a', 'b', 'c'])
        self.assertTrue(self.dm_3x3 == eq_dm)
        self.assertTrue(eq_dm == self.dm_3x3)

    def test_validate(self):
        """Empty stub: SymmetricDistanceMatrix._validate tested elsewhere."""
        pass


class RandomDistanceMatrixTests(TestCase):
    """Tests for bipy.core.distance.random_distance_matrix."""

    def test_default_usage(self):
        """Test generating random distance matrices."""
        exp = DistanceMatrix(np.asarray([[0.0]]), ['1'])
        obs = random_distance_matrix(1)
        self.assertEqual(obs, exp)

        obs = random_distance_matrix(2)
        self.assertEqual(obs.shape, (2, 2))
        self.assertEqual(obs.ids, ('1', '2'))

        obs1 = random_distance_matrix(5)
        num_trials = 10
        found_diff = False
        for _ in range(num_trials):
            obs2 = random_distance_matrix(5)

            if obs1 != obs2:
                found_diff = True
                break

        self.assertTrue(found_diff)

    def test_ids(self):
        """Test generating random dist mats with specific IDs."""
        ids = ['foo', 'bar', 'baz']
        obs = random_distance_matrix(3, ids=ids)
        self.assertEqual(obs.shape, (3, 3))
        self.assertEqual(obs.ids, tuple(ids))

    def test_constructor(self):
        """Test generating random dist mats with a specific constructor."""
        exp = SymmetricDistanceMatrix(np.asarray([[0.0]]), ['1'])
        obs = random_distance_matrix(1, constructor=SymmetricDistanceMatrix)
        self.assertEqual(obs, exp)
        self.assertTrue(isinstance(obs, SymmetricDistanceMatrix))

    def test_random_fn(self):
        """Test passing a different random function than the default."""
        def myrand(num_rows, num_cols):
            # One dm to rule them all...
            data = np.empty((num_rows, num_cols))
            data.fill(42)
            return data

        exp = DistanceMatrix(np.asarray([[0, 42, 42], [42, 0, 42],
                                         [42, 42, 0]]), ['1', '2', '3'])
        obs = random_distance_matrix(3, random_fn=myrand)
        self.assertEqual(obs, exp)

    def test_invalid_input(self):
        """Test error-handling upon invalid input."""
        # Invalid dimensions.
        with self.assertRaises(DistanceMatrixError):
            _ = random_distance_matrix(0)

        # Invalid dimensions.
        with self.assertRaises(ValueError):
            _ = random_distance_matrix(-1)

        # Invalid number of IDs.
        with self.assertRaises(DistanceMatrixError):
            _ = random_distance_matrix(2, ids=['foo'])


# 1x1:
#     0.0
DM_1x1_F = "\ta\na\t0.0\n"

# 2x2:
#       0.0  0.123
#     0.123    0.0
DM_2x2_F = "\ta\tb\na\t0.0\t0.123\nb\t0.123\t0.0\n"

# 2x2, asymmetric:
#     0.0 1.0
#    -2.0 0.0
DM_2x2_ASYM_F = '\ta\tb\na\t0.0\t1.0\nb\t-2.0\t0.0\n'

# 3x3:
#      0.0   0.01   4.2
#     0.01    0.0  12.0
#      4.2   12.0   0.0
DM_3x3_F = ("\ta\tb\tc\na\t0.0\t0.01\t4.2\nb\t0.01\t0.0\t12.0\n"
            "c\t4.2\t12.0\t0.0\n")

# Extra whitespace-only lines throughout. Also has comments before the header.
DM_3x3_WHITESPACE_F = ['# foo',
                       '      \t \t ',
                       ' #bar',
                       '',
                       '',
                       '\ta\t b \tc',
                       'a  \t0.0\t0.01\t4.2',
                       '     \t',
                       'b\t0.01\t0.0\t12.0',
                       '',
                       '\t     \t',
                       '',
                       'c\t4.2\t12.0\t0.0',
                       '',
                       '   \t ',
                       '\t\t\t',
                       ' ']

# missing data
BAD_DM_F1 = 'a\tb\na\t0\t1\nb\t1'

# mismatched IDs
BAD_DM_F2 = '\ta\tb\nb\t0\t1\na\t1\t0'

# extra data lines
BAD_DM_F3 = '\ta\tb\na\t0\t1\nb\t1\t0\n  \nfoo\n\n\n'

# missing data lines
BAD_DM_F4 = '\ta\tb\na\t0\t1\n  \n'

# no data lines
BAD_DM_F5 = '\ta\tb\n'

# non-hollow
BAD_DM_F6 = '\ta\tb\na\t0\t1\nb\t1\t0.1\n'

if __name__ == '__main__':
    main()
