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
from scipy.spatial.distance import euclidean


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

    def test_str(self):
        self.assertRaises(NotImplementedError,
                          Embedding(self.emb, self.seq).__str__)

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
        rstr = repr(p_emb)
        self.assertTrue('SequenceEmbedding' in rstr)
        self.assertTrue(
            '62' in rstr
        )
        self.assertTrue(
            '10' in rstr
        )
        self.assertTrue('IGKEEIQQRL' in rstr)

    def test_str(self):
        emb, s = self.emb, self.seq
        p_emb = SequenceEmbedding(emb, s)

        self.assertEqual(p_emb.__str__(), s)
        self.assertEqual(p_emb.sequence, s)
        self.assertEqual(str(p_emb.ids.tobytes().decode('ascii')), s)

    def test_bytes(self):
        emb, s = self.emb, self.seq
        p_emb = SequenceEmbedding(emb, s)
        res = p_emb.bytes()
        res_str = str(res.tobytes().decode("ascii"))
        self.assertEqual(res_str, s)

    def test_embedding(self):
        emb, s = self.emb, self.seq
        p_emb = SequenceEmbedding(emb, s)
        self.assertEqual(p_emb.embedding.shape, (62, 10))

    def test_assert_length(self):
        with self.assertRaises(ValueError):
            SequenceEmbedding(self.emb, self.seq + "A")


class SequenceVectorTests(TestCase):

    def setUp(self):
        # Create some sample SequenceVector objects for testing
        self.vector1 = np.array([1, 2, 3])
        self.vector2 = np.array([4, 5, 6])
        self.vector3 = np.array([7, 8, 9])
        self.bad_vector = np.array([7, 8])
        self.sequence_vectors = [SequenceVector(self.vector1, "ACGT"),
                                 SequenceVector(self.vector2, "GCTA"),
                                 SequenceVector(self.vector3, "TTAG")]

    def test_vector(self):
        # Test if the vector attribute is set correctly
        for i, vector in enumerate([self.vector1, self.vector2, self.vector3]):
            self.assertTrue(np.array_equal(self.sequence_vectors[i].vector, vector))

    def test_sequence(self):
        # Test if the sequence attribute is set correctly
        for i, sequence in enumerate(["ACGT", "GCTA", "TTAG"]):
            self.assertEqual(self.sequence_vectors[i].sequence, sequence)

    def test_repr(self):
        # Test if the __repr__ method returns the correct string
        for i, sequence_vector in enumerate(self.sequence_vectors):
            self.assertTrue(
                sequence_vector.__repr__().startswith("SequenceVector")
            )
            self.assertTrue(
                'vector' in sequence_vector.__repr__()
            )

            # check latent dimension
            self.assertTrue(
                '4' in sequence_vector.__repr__()
            )

    def test_str(self):
        # Test if the __str__ method returns the correct string
        for i, sequence_vector in enumerate(self.sequence_vectors):
            self.assertEqual(str(sequence_vector), sequence_vector.sequence)

    def test_to_numpy(self):
        # Test if to_numpy returns the correct numpy array
        expected_result = np.array([self.vector1, self.vector2, self.vector3])
        result = SequenceVector.to_numpy(self.sequence_vectors)
        self.assertTrue(np.array_equal(result, expected_result))

    def test_to_numpy_raises(self):
        # assert lengths are equal
        arr = [SequenceVector(self.vector1, "ACGT"),
               SequenceVector(self.vector2, "GCTA"),
               SequenceVector(self.bad_vector, "TTAG")]

        with self.assertRaises(ValueError):
            result = SequenceVector.to_numpy(arr)

    def test_to_distance_matrix(self):
        # Test if to_distance_matrix returns a DistanceMatrix object
        distance_matrix = SequenceVector.to_distance_matrix(self.sequence_vectors)
        self.assertEqual(distance_matrix.shape, (3, 3))
        self.assertTrue(all(isinstance(d, float) for d in distance_matrix.condensed_form()))

        d12 = euclidean(self.vector1, self.vector2)
        d13 = euclidean(self.vector1, self.vector3)
        d23 = euclidean(self.vector2, self.vector3)
        exp_dm = DistanceMatrix([[0, d12, d13],
                                 [d12, 0, d23],
                                 [d13, d23, 0]], ids=["ACGT", "GCTA", "TTAG"])
        self.assertTrue(np.allclose(distance_matrix.data, exp_dm.data))
        self.assertEqual(distance_matrix.ids, exp_dm.ids)

    def test_to_dataframe(self):
        # Test if to_dataframe returns a pandas DataFrame object
        dataframe = SequenceVector.to_dataframe(self.sequence_vectors)
        self.assertIsInstance(dataframe, pd.DataFrame)
        self.assertEqual(dataframe.shape, (3, 3))

        exp_df = pd.DataFrame([self.vector1, self.vector2, self.vector3],
                              index=["ACGT", "GCTA", "TTAG"])
        pd.testing.assert_frame_equal(dataframe, exp_df)

    def test_to_ordination(self):
        # Test if to_ordination returns an OrdinationResults object
        ordination_results = SequenceVector.to_ordination(self.sequence_vectors)
        self.assertEqual(ordination_results.samples.shape, (3, 3))
        self.assertEqual(ordination_results.features.shape, (3, 3))


if __name__ == '__main__':
    main()
