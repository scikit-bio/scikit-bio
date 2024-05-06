# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import copy
import io
import string
from unittest import TestCase, main
from functools import partial
from pathlib import Path
from skbio.util import get_data_path
from skbio.embedding._embedding import (
    SequenceVector, SequenceEmbedding, Embedding)
import numpy as np
import numpy.testing as npt
import pandas as pd
from skbio import DistanceMatrix, OrdinationResults


class EmbeddingTests(TestCase):
    def setUp(self):
        self.emb = np.random.randn(62, 10)
        self.seq = ("IGKEEIQQRLAQFVDHWKELKQLAAARGQRLEESLEYQ"
                    "QFVANVEEEEAWINEKMTLVASED")

    def test_id(self):
        emb, s = self.emb, self.seq
        p_emb = Embedding(emb, list(s))
        np.testing.assert_array_equal(p_emb.ids, np.array(list(s)))

    def test_embedding(self):
        emb, s = self.emb, self.seq
        p_emb = Embedding(emb, s)
        self.assertEqual(p_emb.embedding.shape, (62, 10))

    def test_assert_length(self):
        with self.assertRaises(ValueError):
            Embedding(self.emb, self.seq + "A")


class SequenceEmbeddingTests(TestCase):
    def setUp(self):
        self.emb = np.random.randn(62, 10)
        self.seq = ("IGKEEIQQRLAQFVDHWKELKQLAAARGQRLEESLEYQ"
                    "QFVANVEEEEAWINEKMTLVASED")

    def test_repr(self):
        emb, s = self.emb, self.seq
        p_emb = SequenceEmbedding(emb, s)
        self.assertTrue('SequenceEmbedding' in p_emb.__repr__())


    def test_str(self):
        emb, s = self.emb, self.seq
        p_emb = SequenceEmbedding(emb, s)

        self.assertEqual(p_emb.__str__(), s)
        self.assertEqual(p_emb.sequence, s)
        self.assertEqual(str(p_emb.ids.tobytes().decode('ascii')), s)

    def test_embedding(self):
        emb, s = self.emb, self.seq
        p_emb = SequenceEmbedding(emb, s)
        self.assertEqual(p_emb.embedding.shape, (62, 10))

    def test_assert_length(self):
        with self.assertRaises(ValueError):
            SequenceEmbedding(self.emb, self.seq + "A")


class TestSequenceVectorMethods(TestCase):

    def setUp(self):
        # Create some sample SequenceVector objects for testing
        self.vector1 = np.array([1, 2, 3])
        self.vector2 = np.array([4, 5, 6])
        self.vector3 = np.array([7, 8, 9])
        self.sequence_vectors = [SequenceVector(self.vector1, "ACGT"),
                                 SequenceVector(self.vector2, "GCTA"),
                                 SequenceVector(self.vector3, "TTAG")]

    def test_to_numpy(self):
        # Test if to_numpy returns the correct numpy array
        expected_result = np.array([self.vector1, self.vector2, self.vector3])
        result = SequenceVector.to_numpy(self.sequence_vectors)
        self.assertTrue(np.array_equal(result, expected_result))

    def test_to_distance_matrix(self):
        # Test if to_distance_matrix returns a DistanceMatrix object
        distance_matrix = SequenceVector.to_distance_matrix(self.sequence_vectors)
        self.assertEqual(distance_matrix.shape, (3, 3))
        self.assertTrue(all(isinstance(d, float) for d in distance_matrix.condensed_form()))

    def test_to_dataframe(self):
        # Test if to_dataframe returns a pandas DataFrame object
        dataframe = SequenceVector.to_dataframe(self.sequence_vectors)
        self.assertIsInstance(dataframe, pd.DataFrame)
        self.assertEqual(dataframe.shape, (3, 3))

    def test_to_ordination(self):
        # Test if to_ordination returns an OrdinationResults object
        ordination_results = SequenceVector.to_ordination(self.sequence_vectors)
        self.assertEqual(ordination_results.samples.shape, (3, 3))
        self.assertEqual(ordination_results.features.shape, (3, 3))


if __name__ == '__main__':
    main()
