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

from bipy.core.distance import DistanceMatrix
from bipy.util.unit_test import TestCase, main

class DistanceMatrixTests(TestCase):
    """Tests for the DistanceMatrix class."""

    def setUp(self):
        """Set up test data for use in distance matrix unit tests."""
        self.dm_1x1_f = StringIO(DM_1x1_F)
        self.dm_1x1 = DistanceMatrix([[0.0]], ['a'])

        self.dm_2x2_f = StringIO(DM_2x2_F)
        self.dm_2x2 = DistanceMatrix([[0.0, 0.123], [0.123, 0.0]], ['a', 'b'])

        self.dm_strs = [DM_1x1_F, DM_2x2_F]
        self.dm_fs = [self.dm_1x1_f, self.dm_2x2_f]
        self.dms = [self.dm_1x1, self.dm_2x2]
        self.dm_shapes = [(1, 1), (2, 2)]
        self.dm_sizes = [1, 4]
        self.dm_condensed_data = [np.array([]), np.array([0.123])]

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

    def test_to_file(self):
        """Should serialize a DistanceMatrix to file."""
        for dm_str, dm in izip(self.dm_strs, self.dms):
            for memory_efficient in True, False:
                obs_f = StringIO()
                dm.to_file(obs_f, memory_efficient=memory_efficient)
                obs = obs_f.getvalue()
                obs_f.close()

                self.assertEqual(obs, dm_str)

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

    def test_condensed_data(self):
        """Test retrieving the data matrix in condensed form."""
        for dm, condensed in izip(self.dms, self.dm_condensed_data):
            obs = dm.condensed_data()
            self.assertTrue(np.array_equal(obs, condensed))

# 1x1:
#     0.0
DM_1x1_F = "\ta\na\t0.0\n"

# 2x2:
#     0.0   0.123
#     0.123 0.0
DM_2x2_F = "\ta\tb\na\t0.0\t0.123\nb\t0.123\t0.0\n"

if __name__ == '__main__':
    main()
