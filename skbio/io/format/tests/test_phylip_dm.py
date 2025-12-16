# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import io
import tempfile

import numpy as np

from skbio.util import get_data_path
from skbio.io.format.phylip_dm import (
    _matrix_to_phylip,
    _phylip_dm_sniffer,
    _phylip_dm_to_distance_matrix,
    _parse_phylip_dm_raw,
)
from skbio.io import PhylipFormatError
from skbio.stats.distance import DistanceMatrix


class TestSniffer(unittest.TestCase):
    def setUp(self):
        self.positives = [
            get_data_path(e)
            for e in [
                "phylip_dm_valid_96_lt_phylip_amazon.dist",
                "phylip_dm_valid_98_lt_phylip_amazon.dist",
                "phylip_dm_valid_98_sq_phylip_amazon.dist",
            ]
        ]

        self.negatives = [
            get_data_path(e)
            for e in [
                "empty",
                "whitespace_only",
                "phylip_dm_invalid_empty_line_after_header.dist",
                "phylip_dm_invalid_empty_line_before_header.dist",
                "phylip_dm_sq_invalid_empty_line_after_header.dist",
                "phylip_dm_sq_invalid_empty_line_before_header.dist",
                "phylip_dm_invalid_header_too_long.dist",
                "phylip_dm_invalid_no_header.dist",
                "phylip_dm_invalid_wrong_number_dists.dist",
                "phylip_dm_invalid_wrong_number_dists_sq.dist",
                "phylip_dm_invalid_no_dists.dist",
                "phylip_dm_invalid_zero_header.dist",
            ]
        ]

    def test_positives(self):
        for fp in self.positives:
            self.assertEqual(_phylip_dm_sniffer(fp), (True, {}))

    def test_negatives(self):
        for fp in self.negatives:
            self.assertEqual(_phylip_dm_sniffer(fp), (False, {}))


class TestReaders(unittest.TestCase):
    def setUp(self):
        self.expected_square_matrix = np.loadtxt(get_data_path("dm_raw_data.txt"))

        self.positive_fps = [
            get_data_path(e)
            for e in [
                "phylip_dm_valid_98_lt_phylip_amazon.dist",
                "phylip_dm_valid_98_sq_phylip_amazon.dist",
            ]
        ]

    def test_phylip_dm_to_distance_matrix_valid_files(self):
        for fp in self.positive_fps:
            dm = _phylip_dm_to_distance_matrix(fp)
            self.assertTrue((dm.data == self.expected_square_matrix).all())

    def test_error_empty_file(self):
        with self.assertRaises(PhylipFormatError) as e:
            with open(get_data_path("empty"), "r") as f:
                data = _parse_phylip_dm_raw(f)
        self.assertEqual(str(e.exception), "This file is empty.")

    def test_error_too_many_seqs(self):
        with self.assertRaises(PhylipFormatError) as e:
            with open(get_data_path("phylip_dm_invalid_too_many_seqs.dist"), "r") as f:
                data = _parse_phylip_dm_raw(f)
        self.assertEqual(
            str(e.exception),
            "The number of sequences is not 8 as specified in the header.",
        )


class TestWriters(unittest.TestCase):
    def setUp(self):
        simple_dm = DistanceMatrix(
            [[0, 1, 2, 3], [1, 0, 4, 5], [2, 4, 0, 6], [3, 5, 6, 0]]
        )
        simple_dm_condensed = DistanceMatrix([1, 2, 3, 4, 5, 6], condensed=True)

        self.objs = [simple_dm, simple_dm_condensed]
        self.fps = [
            get_data_path("phylip_dm_simple_sq.dist"),
            get_data_path("phylip_dm_simple_lt.dist"),
        ]

    def test_write(self):
        for fp, obj in zip(self.fps, self.objs):
            fh = io.StringIO()
            _matrix_to_phylip(obj, fh, "\t")
            obs = fh.getvalue()
            fh.close()
            with io.open(fp) as fh:
                exp = fh.read()

            self.assertEqual(obs, exp)

    def test_decorator_write(self):
        for fp, obj in zip(self.fps, self.objs):
            with tempfile.TemporaryFile(mode="w+t") as fh:
                obj.write(fh, format="phylip_dm")
                fh.seek(0)
                obs = fh.read()
            with io.open(fp) as fh:
                exp = fh.read()

            self.assertEqual(obs, exp)

    def test_not_enough_sequences(self):
        with self.assertRaises(PhylipFormatError) as e:
            fh = io.StringIO()
            _matrix_to_phylip(DistanceMatrix([]), fh, "\t")
            fh.close()
        self.assertEqual(
            str(e.exception),
            "DistanceMatrix can only be written in PHYLIP format if there are at least"
            " two samples in the matrix.",
        )


if __name__ == "__main__":
    unittest.main()
