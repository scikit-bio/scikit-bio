# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt

from skbio.alignment import TabularMSA
from skbio.sequence import Sequence, DNA, Protein
from skbio.sequence.distance import p_dist
from skbio.stats.distance import DistanceMatrix

from skbio.alignment._distance import align_dists, _valid_hash


class AlignDistsTests(TestCase):
    def setUp(self):
        self.msa1 = TabularMSA([
            DNA("ATC-GTATCGY"),  # degenerate char: "Y"
            DNA("ATGCG--CCGC"),
            DNA("GT-CGTACG-T"),
            DNA("GT-NTTACAGT"),  # degenerate char: "N"
        ], index=list("abcd"))

        self.msa2 = TabularMSA([
            Protein("MQVEGATSI"),
            Protein("MKNJ-PTSL"),  # degenerate char: "J"
            Protein("MKDE-PYRX"),  # degenerate char: "X"
            Protein("MQUENPT--"),  # non-canonical char: "U"
        ], index=list("abcd"))

    def test_align_dists(self):
        # p-distance, complete deletion (5 sites)
        obs = align_dists(self.msa1, "p_dist")
        self.assertIsInstance(obs, DistanceMatrix)
        self.assertTupleEqual(obs.ids, tuple("abcd"))
        exp = np.array([[0. , 0.2, 0.6, 0.8],
                        [0.2, 0. , 0.4, 0.6],
                        [0.6, 0.4, 0. , 0.4],
                        [0.8, 0.6, 0.4, 0. ]])
        npt.assert_array_equal(obs.data, exp)

        # supply function
        obs = align_dists(self.msa1, p_dist)
        npt.assert_array_equal(obs.data, exp)

        # pairwise deletion
        obs = align_dists(self.msa1, "p_dist", shared_by_all=False)
        exp = np.array([[0. , 2/7, 3/7, 4/8],
                        [2/7, 0. , 3/7, 4/7],
                        [3/7, 3/7, 0. , 2/8],
                        [4/8, 4/7, 2/8, 0. ]])
        npt.assert_allclose(obs.data, exp)

        # supply function
        obs = align_dists(self.msa1, p_dist, shared_by_all=False)
        npt.assert_array_equal(obs.data, exp)


    def test_valid_hash(self):
        # DNA sequences
        seqs = np.vstack([seq._bytes for seq in self.msa1])

        obs = _valid_hash("ACGT", DNA)[seqs]
        exp = np.array([[1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0],
                        [1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1],
                        [1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1],
                        [1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1]], dtype=bool)
        npt.assert_array_equal(obs, exp)

        obs = _valid_hash("ACGTN", DNA)[seqs]
        exp = np.array([[1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0],
                        [1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1],
                        [1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1],
                        [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1]], dtype=bool)
        npt.assert_array_equal(obs, exp)

        obs = _valid_hash("ABCD", DNA)[seqs]
        exp = np.array([[1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0],
                        [1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1],
                        [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0]], dtype=bool)
        npt.assert_array_equal(obs, exp)

        obs = _valid_hash("nongap", DNA)[seqs]
        exp = np.array([[1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1],
                        [1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1],
                        [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1]], dtype=bool)
        npt.assert_array_equal(obs, exp)

        obs = _valid_hash("definite", DNA)[seqs]
        exp = np.array([[1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0],
                        [1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1],
                        [1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1],
                        [1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1]], dtype=bool)
        npt.assert_array_equal(obs, exp)

        # protein sequences
        seqs = np.vstack([seq._bytes for seq in self.msa2])

        obs = _valid_hash("nongap", Protein)[seqs]
        exp = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 0, 1, 1, 1, 1],
                        [1, 1, 1, 1, 0, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 0, 0]], dtype=bool)
        npt.assert_array_equal(obs, exp)

        obs = _valid_hash("definite", Protein)[seqs]
        exp = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 0, 0, 1, 1, 1, 1],
                        [1, 1, 1, 1, 0, 1, 1, 1, 0],
                        [1, 1, 1, 1, 1, 1, 1, 0, 0]], dtype=bool)
        npt.assert_array_equal(obs, exp)

        obs = _valid_hash("canonical", Protein)[seqs]
        exp = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 0, 0, 1, 1, 1, 1],
                        [1, 1, 1, 1, 0, 1, 1, 1, 0],
                        [1, 1, 0, 1, 1, 1, 1, 0, 0]], dtype=bool)
        npt.assert_array_equal(obs, exp)


if __name__ == "__main__":
    main()
