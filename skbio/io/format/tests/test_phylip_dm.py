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
                "phylip_dm_valid_lt.dist",
                "phylip_dm_valid_sq.dist",
                "phylip_dm_good_simple_strict_square.dist",
                "phylip_dm_simple_lt.dist",
                "phylip_dm_simple_sq.dist"
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
                "phylip_dm_invalid_no_dists.dist",
                "phylip_dm_invalid_zero_header.dist",
                "phylip_dm_invalid_too_many_columns.dist",
                "phylip_dm_invalid_too_few_columns_sq.dist",
                "phylip_dm_invalid_wrong_number_dists_lt.dist",
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
        self.expected_data_relaxed = [
            np.array(
                [
                    [0.0, 1.7043, 2.0235, 2.1378, 1.5232, 1.8261, 1.9182, 2.0039, 1.9431, 1.9663, 2.0593, 1.6664, 1.732, 1.7101,],
                    [1.7043, 0.0, 1.1901, 1.3287, 1.2423, 1.2508, 1.2536, 1.3066, 1.2827, 1.3296, 1.2005, 1.346, 1.3757, 1.3956,],
                    [2.0235, 1.1901, 0.0, 1.2905, 1.3199, 1.3887, 1.4658, 1.4826, 1.4502, 1.8708, 1.5356, 1.4577, 1.7803, 1.6661,],
                    [2.1378, 1.3287, 1.2905, 0.0, 1.7878, 1.3137, 1.3788, 1.3826, 1.4543, 1.6683, 1.6606, 1.5935, 1.7119, 1.7599,],
                    [1.5232, 1.2423, 1.3199, 1.7878, 0.0, 1.0642, 1.1124, 0.9832, 1.0629, 0.9228, 1.0681, 0.9127, 1.0635, 1.0557,],
                    [1.8261, 1.2508, 1.3887, 1.3137, 1.0642, 0.0, 0.1022, 0.2061, 0.3895, 0.8035, 0.7239, 0.7278, 0.7899, 0.6933,],
                    [1.9182, 1.2536, 1.4658, 1.3788, 1.1124, 0.1022, 0.0, 0.2681, 0.393, 0.7109, 0.729, 0.7412, 0.8742, 0.7118,],
                    [2.0039, 1.3066, 1.4826, 1.3826, 0.9832, 0.2061, 0.2681, 0.0, 0.3665, 0.8132, 0.7894, 0.8763, 0.8868, 0.7589,],
                    [1.9431, 1.2827, 1.4502, 1.4543, 1.0629, 0.3895, 0.393, 0.3665, 0.0, 0.7858, 0.714, 0.7966, 0.8288, 0.8542,],
                    [1.9663, 1.3296, 1.8708, 1.6683, 0.9228, 0.8035, 0.7109, 0.8132, 0.7858, 0.0, 0.7095, 0.5959, 0.6213, 0.5612,],
                    [2.0593, 1.2005, 1.5356, 1.6606, 1.0681, 0.7239, 0.729, 0.7894, 0.714, 0.7095, 0.0, 0.4604, 0.5065, 0.47,],
                    [1.6664, 1.346, 1.4577, 1.5935, 0.9127, 0.7278, 0.7412, 0.8763, 0.7966, 0.5959, 0.4604, 0.0, 0.3502, 0.3097,],
                    [1.732, 1.3757, 1.7803, 1.7119, 1.0635, 0.7899, 0.8742, 0.8868, 0.8288, 0.6213, 0.5065, 0.3502, 0.0, 0.2712,],
                    [1.7101, 1.3956, 1.6661, 1.7599, 1.0557, 0.6933, 0.7118, 0.7589, 0.8542, 0.5612, 0.47, 0.3097, 0.2712, 0.0,],
                ]
            ),
            np.array(
                [
                    [0.0, 1.0, 2.0, 3.0, 3.0],
                    [1.0, 0.0, 2.0, 3.0, 3.0],
                    [2.0, 2.0, 0.0, 3.0, 3.0],
                    [3.0, 3.0, 3.0, 0.0, 1.0],
                    [3.0, 3.0, 3.0, 1.0, 0.0],
                ]
            ),
        ]

        self.positive_fps_relaxed = [
            get_data_path(e)
            for e in [
                "phylip_dm_valid_lt.dist",
                "phylip_dm_valid_sq.dist",
            ]
        ]

    def test_phylip_dm_to_distance_matrix_valid_files(self):
        for fp, exp in zip(self.positive_fps_relaxed, self.expected_data_relaxed):
            dm = _phylip_dm_to_distance_matrix(fp)
            self.assertTrue((dm.data == exp).all())

    def test_error_empty_file(self):
        with self.assertRaises(PhylipFormatError) as e:
            with open(get_data_path("empty"), "r") as f:
                data = _parse_phylip_dm_raw(f)
        self.assertEqual(str(e.exception), "This file is empty.")

    def test_error_wrong_number_seqs(self):
        fps = ["phylip_dm_invalid_too_many_rows.dist",
               "phylip_dm_invalid_too_few_rows.dist"]
        for fp in fps:
            with self.assertRaises(PhylipFormatError) as e:
                with open(get_data_path(fp), "r") as f:
                    data = _parse_phylip_dm_raw(f)
            self.assertEqual(
                str(e.exception),
                "The number of sequences is not 5 as specified in the header.",
            )

    def test_error_matrix_data_parsed_as_id(self):
        fp = "phylip_dm_good_simple_strict_square_reader.dist"
        with self.assertRaises(PhylipFormatError) as e:
            with open(get_data_path(fp), "r") as f:
                dm = DistanceMatrix.read(f, format="phylip_dm", strict=False)
        self.assertEqual(
            str(e.exception),
            "Inconsistent distance counts detected: [4, 4, 4, 3]. This may indicate that sequence IDs contain whitespace. IDs may only contain whitespace if the strict parameter is set to True. Expected either all 4 (square) or 0,1,2,... (lower triangular)."
            )


class TestWriters(unittest.TestCase):
    def setUp(self):
        self.dm_condensed = DistanceMatrix([1, 2, 3, 4, 5, 6], condensed=True)
        self.dm_square = DistanceMatrix(
            [[0, 1, 2, 3], [1, 0, 4, 5], [2, 4, 0, 6], [3, 5, 6, 0]]
        )
        
        self.expected_lower_tri = "4\n0\n1\t1.0\n2\t2.0\t4.0\n3\t3.0\t5.0\t6.0\n"
        self.expected_square = "4\n0\t0.0\t1.0\t2.0\t3.0\n1\t1.0\t0.0\t4.0\t5.0\n2\t2.0\t4.0\t0.0\t6.0\n3\t3.0\t5.0\t6.0\t0.0\n"
    
    def test_condensed_to_lower_tri(self):
        """Test writing condensed DistanceMatrix to lower triangular format."""
        fh = io.StringIO()
        self.dm_condensed.write(fh, format='phylip_dm', lower_tri=True)
        result = fh.getvalue()
        fh.close()
        
        self.assertEqual(result, self.expected_lower_tri)
    
    def test_condensed_to_square(self):
        """Test writing condensed DistanceMatrix to square format."""
        fh = io.StringIO()
        self.dm_condensed.write(fh, format='phylip_dm', lower_tri=False)
        result = fh.getvalue()
        fh.close()
        
        self.assertEqual(result, self.expected_square)
    
    def test_square_to_lower_tri(self):
        """Test writing square DistanceMatrix to lower triangular format."""
        fh = io.StringIO()
        self.dm_square.write(fh, format='phylip_dm', lower_tri=True)
        result = fh.getvalue()
        fh.close()
        
        self.assertEqual(result, self.expected_lower_tri)
    
    def test_square_to_square(self):
        """Test writing square DistanceMatrix to square format."""
        fh = io.StringIO()
        self.dm_square.write(fh, format='phylip_dm', lower_tri=False)
        result = fh.getvalue()
        fh.close()
        
        self.assertEqual(result, self.expected_square)

    def test_not_enough_sequences(self):
        with self.assertRaises(PhylipFormatError) as e:
            fh = io.StringIO()
            _matrix_to_phylip(DistanceMatrix([]), fh, "\t", lower_tri=True)
            fh.close()
        self.assertEqual(
            str(e.exception),
            "DistanceMatrix can only be written in PHYLIP format if there are at least"
            " two samples in the matrix.",
        )


if __name__ == "__main__":
    unittest.main(buffer=False)
