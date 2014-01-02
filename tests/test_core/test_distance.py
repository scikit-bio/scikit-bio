#!/usr/bin/env python
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from StringIO import StringIO

from bipy.core.distance import DistanceMatrix
from bipy.util.unit_test import TestCase, main

class DistanceMatrixTests(TestCase):
    """Tests for the DistanceMatrix class."""

    def setUp(self):
        """Set up test data for use in distance matrix unit tests."""
        self.dm_1x1_f = StringIO(DM_1x1_F)
        self.dm_1x1 = DistanceMatrix([[0.0]], ['a'])

    def test_round_trip_read_write(self):
        """Test reading, writing, and reading again works as expected."""
        # Read.
        dm1 = DistanceMatrix.from_file(self.dm_1x1_f)

        # Write.
        out_f = StringIO()
        dm1.to_file(out_f)
        out_f.seek(0)

        # Read.
        dm2 = DistanceMatrix.from_file(out_f)
        self.assertEqual(dm1, dm2)

    def test_from_file(self):
        """Should parse and return a valid DistanceMatrix given a file."""
        obs = DistanceMatrix.from_file(self.dm_1x1_f)
        self.assertEqual(obs, self.dm_1x1)

    def test_to_file(self):
        """Should serialize a DistanceMatrix to file."""
        obs_f = StringIO()
        self.dm_1x1.to_file(obs_f)
        obs = obs_f.getvalue()
        obs_f.close()

        self.assertEqual(obs, DM_1x1_F)

DM_1x1_F = "\ta\na\t0.0\n"

if __name__ == '__main__':
    main()
